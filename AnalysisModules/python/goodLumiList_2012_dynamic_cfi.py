import FWCore.ParameterSet.Config as cms

#--- This configuration uses the JSON file:
#--- /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-191859_8TeV_PromptReco_Collisions12_JSON.txt
lumisToProcess = cms.untracked.VLuminosityBlockRange()
lumisToProcess.extend([
   '190645:10-190645:110',
   '190704:1-190704:3',
   '190705:1-190705:5',
   '190705:7-190705:76',
   '190705:78-190705:184',
   '190738:1-190738:130',
   '190738:133-190738:226',
   '190738:229-190738:355',
   '191057:1-191057:1',
   '191057:4-191057:65',
   '191062:1-191062:1',
   '191062:3-191062:3',
   '191062:5-191062:214',
   '191062:216-191062:549',
   '191090:1-191090:331',
   '191201:38-191201:49',
   '191201:52-191201:79',
   '191202:1-191202:64',
   '191202:66-191202:68',
   '191202:87-191202:105',
   '191202:108-191202:118',
   '191226:77-191226:78',
   '191226:81-191226:831',
   '191226:833-191226:1454',
   '191226:1456-191226:1466',
   '191226:1469-191226:1507',
   '191226:1510-191226:1686',
   '191247:1-191247:153',
   '191247:156-191247:280',
   '191247:283-191247:606',
   '191247:608-191247:620',
   '191247:622-191247:818',
   '191247:821-191247:834',
   '191247:837-191247:1031',
   '191247:1034-191247:1046',
   '191247:1049-191247:1140',
   '191247:1143-191247:1187',
   '191247:1190-191247:1214',
   '191247:1217-191247:1224',
   '191248:1-191248:103',
   '191264:59-191264:79',
   '191264:82-191264:152',
   '191264:155-191264:189',
   '191271:56-191271:158',
   '191276:1-191276:16',
   '191277:1-191277:28',
   '191277:30-191277:164',
   '191277:167-191277:253',
   '191277:255-191277:457',
   '191277:460-191277:535',
   '191277:537-191277:576',
   '191277:579-191277:775',
   '191277:778-191277:811',
   '191277:813-191277:849',
   '191367:1-191367:276',
   '191411:1-191411:23',
   '191695:1-191695:1',
   '191697:3-191697:9',
   '191718:43-191718:95',
   '191718:98-191718:207',
   '191720:1-191720:1',
   '191720:3-191720:15',
   '191720:17-191720:181',
   '191721:1-191721:1',
   '191721:3-191721:34',
   '191721:36-191721:183',
   '191721:186-191721:189',
   '191726:1-191726:13',
   '191810:15-191810:15',
   '191810:22-191810:49',
   '191810:52-191810:92'
])
