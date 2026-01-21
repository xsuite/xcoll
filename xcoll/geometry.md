PSEUDOCODE for material interactions tracking in collimator
===========================================================

- create DriftTrajectory from particle
- loop:
    - find first crossing with collimator shape (segments)
    - if distance larger than end of BeamElement -> move particle to exendit and break loop
    - move particle to collimator impact point
    - create MCSTrajectory from new particle
    - loop:
        - find first crossing with collimator shape (segments)
        - do not move, but calculate path length
        - compare to distance to nuclear interaction
        - shortest wins (nuclear or exit)
        - move to next point
    - when at exit, create DriftTrajectory from particle



PSEUDOCODE for material interactions tracking in crystal
========================================================

- create DriftTrajectory from particle
- loop:
    - find first crossing with crystal shape (segments)
    - if distance larger than end of BeamElement -> move particle to end and break loop
    - move particle to crystal impact point
    - create MCSTrajectory (if amorphous) or CircularTrajectory (if channelling) from new particle
    - loop:
        - find first crossing with crystal shape (segments)
        - do not move, but calculate path length
        - compare to distance to nuclear interaction and dechannelling
        - shortest wins (nuclear, dechannelling, or exit)
        - move to next point
    - when at exit, create DriftTrajectory from particle


PSEUDOCODE for crossing finder
==============================

- if DriftTrajectory and any segment
=> use analytical