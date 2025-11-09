# DD-BOAT 5 : Guerlédan Partie 1

Ce respository présente un système de navigation autonome pour bateau robotique (DD-Boat) utilisant un GPS, une IMU et contrôle feedforward avec correction de cap réalisé dans le cadre de la semaine à Guerlédan de l'ENSTA campus de Brest.

## 📋 Description

Ce projet implémente un système de navigation autonome permettant à un bateau de suivre une trajectoire prédéfinie en temps réel. Le système utilise:
- **GPS** pour la localisation
- **IMU 9 axes** (magnétomètre, accéléromètre, gyroscope) pour l'orientation
- **Contrôle feedforward** avec correction de cap pour le pilotage
- **Calibration magnétomètre** pour compenser les distorsions

Le scénario principal consiste à suivre une "proie" virtuelle se déplaçant sur un cercle, avec le bateau positionné sur une orbite décroissante autour de celle-ci.

## 🎯 Fonctionnalités

- Navigation GPS avec conversion en repère local
- Calibration automatique du magnétomètre (correction d'ellipsoïde)
- Suivi de trajectoire temporelle avec anticipation de vitesse
- Génération de traces GPX
- Sauvegarde automatique des données (même en cas d'interruption)
- Synchronisation temporelle pour départs multiples coordonnés
- Simulation et visualisation des trajectoires

## 📦 Structure du projet

```
.
├── Autonomous_capture.py         # Script principal
├── readIMU.py                    # Acquisition des données IMU pour calibration
├── plotIMUCalib.py              # Visualisation de la calibration magnétomètre
├── Simulateur_Trajectory_Simple.py  # Simulation des trajectoires sans animations
├── drivers-ddboat-v2/           # Drivers matériels
│   ├── arduino_driver_v2.py     # Communication Arduino (moteurs)
│   ├── imu9_driver_v2.py        # Driver IMU 9 axes
│   └── gps_driver_v2.py         # Driver GPS
└── datacalibration.txt          # Données de calibration (générés par readIMU.py)
```
## 📖 Structure du script principal

### Vue d'ensemble de `Autonomous_capture.py`

Le script est organisé en modules fonctionnels hiérarchiques:

#### 1. **Utilitaires mathématiques**

```python
sur_360(a)              # Normalise un angle dans [0, 360°]
ang_err(target, meas)   # Erreur angulaire dans [-180°, 180°]
sawtooth(a)             # Normalise un angle dans [-π, π]
bearing_geo_from_vec()  # Convertit vecteur (vx,vy) en cap géographique
```

#### 2. **Conversion GPS (NMEA → degrés décimaux)**

```python
gll_ddmm_to_dd(st)     # Convertit format DDMM.MMMM en degrés décimaux
local_xy()             # Projette (lat,lon) en coordonnées locales (x,y)
```

Utilise une projection locale tangente avec rayon terrestre R=6378137m.

#### 3. **Calibration magnétomètre**

```python
fit_ellispoide_sphere(M)    # Ajuste ellipsoïde par SVD
load_mag_calibration(path)  # Charge calibration depuis CSV
```

**Algorithme:**
- Résolution par décomposition SVD: D·p = 0
- Extraction des paramètres de l'ellipsoïde (Q, a, c)
- Calcul du centre b et de la matrice de transformation A
- Teste les deux signes de p pour trouver une solution positive définie

#### 4. **Classe `FFCtrl` - Contrôleur Feedforward**

```python
class FFCtrl:
    __init__(imu, ard, A_mag, b_mag, alpha_v, alpha_w, kp_psi, mmax, lam)
    tanh_sat(x)              # Saturation tanh des commandes
    read_heading_rad()        # Lecture cap calibré
    step(u_d, omega_d, psi_d) # Calcul et envoi commandes moteurs
```

**Loi de commande:**
```
ψ = read_heading_rad()
e_ψ = sawtooth(ψ - ψ_d)
ω_cmd = ω_d - kp_psi · e_ψ

left  = tanh_sat(alpha_v·u_d + alpha_w·ω_cmd)
right = tanh_sat(alpha_v·u_d - alpha_w·ω_cmd)
```

#### 5. **Navigation GPS**

```python
class NavState:          # État: (x, y, lat, lon)
poll_nav(gps)            # Lecture non-bloquante GPS
wait_first_fix(gps)      # Attente signal GPS valide
goto_point()             # Navigation vers un point (x,y)
```

**Stratégie `goto_point`:**
1. Calcul distance et direction vers objectif
2. Si distance ≤ rayon → arrivé
3. Sinon: calcul cap désiré ψ_d, application FFC
4. Sauvegarde point GPX à chaque mise à jour

#### 6. **Classe `PreyOrbit` - Générateur de trajectoire**

```python
class PreyOrbit:
    __init__(Cx, Cy, R1, w1, N, kidx)
    R2(t)                    # Rayon orbite décroissant: 10·exp(-t/200) + 5
    w2(t)                    # Vitesse angulaire: 1/(2·R2(t))
    proie(t)                 # Position proie: C + R1·[cos(w1·t), sin(w1·t)]
    desired_pos(t)           # Position désirée bateau
    desired_vel(t, dt)       # Vitesse désirée (dérivée numérique)
```

**Formule trajectoire bateau:**
```
a(t) = p(t) + R2(t)·[cos(θ2), sin(θ2)]
où θ2 = w2·t + 2π·kidx/N
```

#### 7. **Fonction `track_at()` - Suivi de trajectoire**

Boucle principale de guidage:

```python
def track_at(ctrl, gps, traj, t_start, u_d_nom, k_pos, dt_cmd, duration_s):
    while time.time() - t_begin < duration_s:
        # 1. Lecture position GPS
        ns = poll_nav(gps)
        
        # 2. Calcul trajectoire désirée
        t = time.time() - t_start
        xd, yd = traj.desired_pos(t)
        vxd, vyd = traj.desired_vel(t)
        
        # 3. Loi de guidage par anticipation
        ex, ey = xd - ns.x, yd - ns.y
        vrefx = vxd + k_pos · ex
        vrefy = vyd + k_pos · ey
        
        # 4. Calcul cap et vitesse angulaire désirés
        psi_d = bearing_geo_from_vec(vrefx, vrefy)
        omega_d = (psi_d - psi_d_prev) / dt
        
        # 5. Application du contrôle
        ctrl.step(u_d=u_d_nom, omega_d=omega_d, psi_d=psi_d)
```

**Principe:** La loi de guidage anticipe la position future en corrigeant l'erreur de position avec un gain k_pos.

#### 8. **Fonction `main()` - Orchestration**

Séquence d'exécution:

```python
1. Parsing arguments (argparse)
2. Initialisation matériel (Arduino, IMU, GPS)
3. Chargement calibration magnétomètre
4. Création contrôleur FFCtrl
5. Initialisation GPX
6. Calcul paramètres trajectoire (centre local, PreyOrbit)
7. Définition point HOME
8. [OPTIONNEL] Navigation vers point de départ
9. Attente synchronisation temporelle (wait_for)
10. Installation handler SIGINT (Ctrl+C → sauvegarde GPX)
11. SUIVI: track_at() pendant duration_s
12. Arrêt moteurs
13. Sauvegarde trace GPX
```

#### 9. **Gestion des interruptions (SIGINT)**

```python
def _on_sigint(sig, frame):
    # Sauvegarde GPX de manière atomique
    with open(args.gpx, 'w') as f:
        f.write(gpx.to_xml())
        f.flush()
        os.fsync(f.fileno())
    # Arrêt moteurs
    ctrl.ard.send_arduino_cmd_motor(0, 0)
    sys.exit(130)
```

**Sécurité:** Garantit la sauvegarde des données même en cas d'arrêt brutal (Ctrl+C, perte de connexion SSH).

### Flux de données

```
GPS → poll_nav() → NavState(x, y, lat, lon)
                      ↓
IMU → read_heading_rad() → ψ (cap calibré)
                      ↓
PreyOrbit → desired_pos(t), desired_vel(t)
                      ↓
track_at() → Loi de guidage → (u_d, ω_d, ψ_d)
                      ↓
FFCtrl.step() → Commandes moteurs (left, right)
                      ↓
Arduino → Actionneurs
```

### Points techniques importants

**Gestion du temps:**
- `wait_for(hh:mm:ss)`: synchronisation multi-bateaux
- `time.time() - t_start`: temps relatif pour trajectoire

**Robustesse GPS:**
- Lecture non-bloquante avec `poll_nav()`
- Réutilisation dernière position valide si échec lecture
- Timeout et validation distance minimale

**Calibration magnétomètre:**
- Automatique si `datacalibration.txt` existe
- Sinon: utilisation données brutes (moins précis)

**Format GPX:**
- Un track avec un segment par mission
- Ajout point à chaque mise à jour GPS (~0.1s)

### Configuration matérielle

Le système nécessite:
- Arduino pour le contrôle moteur
- IMU 9 axes (magnétomètre, accéléromètre, gyroscope)
- Module GPS
- Drivers dans le dossier `drivers-ddboat-v2/`

## 🚀 Utilisation

### 1. Calibration du magnétomètre

Effectuer une calibration avant la première utilisation:

```bash
# Acquisition des données (tourner le bateau dans toutes les directions)
python readIMU.py

# Visualisation de la calibration
python plotIMUCalib.py
```

Cela génère `datacalibration.txt` utilisé automatiquement par le script principal.

### 2. Simulation des trajectoires

Avant le déploiement, visualiser les trajectoires:

```bash
python Simulateur_Trajectory_Simple.py
```

Ajuster les paramètres dans le fichier selon vos besoins.

### 3. Navigation autonome

#### Utilisation basique

```bash
python Autonomous_capture.py
```

#### Options principales

```bash
# Centre de la trajectoire (latitude, longitude)
--C_latlon 48.199706 -3.018784

# Rayon du cercle de la proie
--R1 240

# Vitesse angulaire de la proie (rad/s)
--w1 -0.0013889

# Nombre de bateaux (pour missions coordonnées)
--N 9

# Index du bateau (détermine sa position sur l'orbite)
--kidx 5

# Commande moteur normalisée
--u_d 2.0

# Gains de contrôle
--kp_psi 1.2
--k_pos 0.8

# Durée de la mission
--duration 750

# Départ synchronisé à une heure précise
--start_time "14:30:00"

# Fichier GPX de sortie
--gpx trace_nav.gpx
```

#### Exemple complet

```bash
python Autonomous_capture.py \
  --C_latlon 48.199706 -3.018784 \
  --R1 240 \
  --N 9 \
  --kidx 5 \
  --u_d 2.0 \
  --duration 750 \
  --start_time "10:00:00" \
  --gpx mission_boat5.gpx
```

## 📐 Configuration du repère

Modifier les constantes dans `Autonomous_capture.py`:

```python
REF_LAT = 48.1991683  # Latitude de référence
REF_LON = -3.01473    # Longitude de référence
```

Ces coordonnées définissent l'origine (0,0) du repère local.

## 🎮 Contrôle

Le système utilise un **contrôleur feedforward avec correction de cap** (FFC):

- **u_d**: commande longitudinale (vitesse avant)
- **ω_d**: vitesse angulaire désirée
- **ψ_d**: cap désiré
- Correction proportionnelle sur l'erreur de cap
- Saturation tanh pour limiter les commandes moteurs

### Paramètres de réglage

- `alpha_v`: coefficient vitesse linéaire (défaut: 1.0)
- `alpha_w`: coefficient vitesse angulaire (défaut: 0.7)
- `kp_psi`: gain de correction de cap (défaut: 1.2)
- `mmax`: limite de commande moteur (défaut: 200)
- `lam`: coefficient de saturation tanh (défaut: 3.0)

## 📊 Modèle de trajectoire

### Trajectoire de la proie

```
p(t) = C + R1 · [cos(w1·t), sin(w1·t)]
```

- **C**: centre du cercle
- **R1**: rayon fixe
- **w1**: vitesse angulaire constante

### Trajectoire du bateau

```
a(t) = p(t) + R2(t) · [cos(θ2), sin(θ2)]
```

Avec:
- **R2(t) = 10·exp(-t/200) + 5**: rayon décroissant
- **w2(t) = 1/(2·R2(t))**: vitesse angulaire adaptative
- **θ2 = φ0 + w2·t**: angle avec déphasage initial φ0 = 2π·k/N

## 📁 Fichiers générés

- `datacalibration.txt`: données de calibration magnétomètre
- `trace_nav.gpx`: trace GPS de la mission
- `nuage_avant_xyz.png`: visualisation magnétomètre avant calibration
- `nuage_apres_xyz.png`: visualisation magnétomètre après calibration

## 🛡️ Sécurité

- **Arrêt d'urgence**: Ctrl+C sauvegarde automatiquement la trace GPX et arrête les moteurs
- **Timeouts**: protection contre la perte de signal GPS
- **Saturation moteurs**: limitation des commandes pour éviter les instabilités

## 🔬 Algorithmes

### Calibration magnétomètre

Ajustement d'ellipsoïde par moindres carrés (SVD) pour corriger:
- Hard iron (offset magnétique)
- Soft iron (déformation du champ)

Transformation: **X_calibré = A · (X_raw - b)**

### Suivi de trajectoire

Loi de guidage par anticipation:

```
v_ref = v_desired + k_pos · (x_desired - x_actual)
ψ_d = atan2(v_ref_y, v_ref_x)
```

## Crédits

Scripts réalisés par BIECHE Matys Adéas et RAMIS Lancelot.
Les parties où ChatGPT a été utilisé sont mis en évidence. 
Utilise les drivers DDBoat v2.