# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import (Material, CrystalMaterial, RefMaterial,
                       _DEFAULT_MATERIAL, _DEFAULT_CRYSTALMATERIAL)
from .atoms import *
from .allotropes import *
from .compounds import *
from .mixtures import *
from .sixtrack import *
from .crystals import *

# Freeze all material instances (to prevent accidental modification)
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and name != 'obj':
        obj._frozen = True
del name, obj

# Get database instance and expose some methods
from .database import db
show = db.show

# Some sources:
# https://epaper.kek.jp/ipac2023/pdf/WEPA148.pdf
# https://pdf.sciencedirectassets.com/271508/1-s2.0-S0008622318X00074/1-s2.0-S0008622318303555/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJD%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIEpatVR2iMAynQtBCSVxWwxx0WFnuPyFK%2B3jgOJ2INnhAiBBw6VGFuLSEoweOC7ryGAImpGk8ApFmp%2BBYwt3hQFXEyqyBQg5EAUaDDA1OTAwMzU0Njg2NSIMcqyr5hKtD3HLV7%2F2Ko8FQOHztlojw87Oi7LV2xR%2Bkwk9JkxSumwlU4aqcgWrDL%2BHnheO7yujorf9v3HMoJE906w4h8ph2cj2KbX7USxL7x0aHkDG8Sm0%2Fx7lRohiq3oQtWF4WucsVEvai%2BiyKff%2FTdUNFhaR%2BICPjMA4Tl%2F5YepW6lgDP1S7UBwCkxl4ymSdbCS8Hbtevx%2FyJIYoEuZHinZfv6yCX0w1iqlvCCaUb2m9haj%2BCHQjCBtqClnQF463e%2FC4j8EApMoidsggxBGXEqNlc%2FuxAcJz4YJuv5a5dK4L9g3mOhBaUFLf5tmSmMptQG9BlIr8MdVCsrOdjbn5u0ho9g%2BWCzEap0YR%2FUKmZ7PkcmZalTPTglyBy%2BouXmiJJ1QSU7mPw%2BXEizno8572zjv8eA6a7aziCwG3tHjKXabe%2FQrHu%2BzoxkF5TgpzyPPuWLxf5XAGrRMBrXju4b8u1eYsQQTAelTBlxtT9WSg9J0M1fFAL%2Bq870p0OtZsTwMWqN6025%2FD6FbwmItW%2FTS4GIBS%2B2E0nRIQvZn2mq3HL3w2okPFjEZVluR1q3vH5NTUB9KSFeQROFjVUoC4hDHKE3fgxX8hfPboKFbN73MrkMI65x0vj42PrqmbQvWHxV16ZtOrBVMKfiIW5ffK6te9Zyi%2FCC7%2BFQKwOfZg6J7O3pSIN5tKQDmXwV%2FuxAvYEK2iLJzh3s6oMSOZdYWeKoy1of3D0TzvajHCThM99hLWt19YEJzj%2B2EE9jUfbRv1n1Ctu%2FKEQpYYYV2qF7Yvw7ivwscTOrDikRjnY%2F2SE6C6kOynmgoXoJ3IiDxN5a%2BZto3B4y6ZbRPacmm96Pr5jJSDksq9n9Gb56kRg%2FyPoACGUkhSsVra8aISegF87wPMaDC897DHBjqyAUxY1XGbgOrxrC2NO1l7QQE0tMrmXr23nsZJ0htjCkO7FmG8Q2PreWfasjJuK0EpwwnDEMdzdvG4GSJMngSBLrFrXs7R%2F1%2Fj8gSd4ZuMjkdY1E22qORRSi9MSaDi1mYgEu6W5S9Kq5u5avLz2ZvrT5pQn0rd2JhZpQpl9lyNp67CYGbX5PzYIpKJ5Yw8msrTH5wxt8tagssiB6xHnpgnsDLo%2BDxXVXlwd5FLI9KNB1MlZro%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20251012T235435Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY7NQSSLCB%2F20251012%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=378e6b9cacb299d5ebc48df11583b1db9a56a9c38ea4eb0088fa25181d79bc5e&hash=ba8572e278d3740d984f493cdcc196226e09d391f2bfbd245654826268e741d3&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0008622318303555&tid=spdf-dc8d6dd3-57b0-41ef-b769-ee6273970373&sid=69e545f78f7a2947339ac3f-be9cff1d6e6fgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=05155858020e5e5003&rr=98da85f3cd01bc81&cc=fr
# https://indico.cern.ch/event/1076808/contributions/5047893/attachments/2511128/4316133/Jaw-materials-justifications_v5_NM.xlsx
# https://edms.cern.ch/ui/#!master/navigator/document?D:101081798:101081798:approvalAndComments
# https://iopscience.iop.org/article/10.1088/1742-6596/2687/8/082037/pdf
# https://indico.cern.ch/event/340703/contributions/802090/attachments/668678/919103/2015_26_03_collmat.pdf


# TCP: CarbonFibreCarbon AC150GPH
# TCPPM: MolybdenumGraphite MG6403Fc
# TCSG: AC150GPH
# TCSP: AC150GPH
# TCSPM: MG6403Fc
# TCT: INERM180
# TCL: COPPER
# TCLA: INERM180
# TDI: BN5000 / Al / Cu
# TCLIA: GraphiteR4550 GRAR4550
# TCLIB: AC150GPH
# TCDQ: CC_1_75 / CC_1_40 / CC_1_75
