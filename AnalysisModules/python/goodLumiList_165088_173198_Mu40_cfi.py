import FWCore.ParameterSet.Config as cms

#--- This configuration was updated to include the following ---#
#--- Run 2011A: prompt V4                                    ---#
#--- Run 2011A: re-reco August 5                             ---#
#--- Run 2011A: prompt V6                                    ---#
#--- Including ONLY lumi sections with HLT_Mu40 unprescaled  ---#
#--- There are more prompt V6 events, but not here...        ---#
lumisToProcess = cms.untracked.VLuminosityBlockRange()
lumisToProcess.extend([
        #--- Prompt V4 ---#
	'160431:19-160431:218',
	'160577:254-160577:306',
	'160578:6-160578:53',
	'160578:274-160578:400',
	'160871:68-160871:208',
	'160872:1-160872:9',
	'160872:25-160872:35',
	'160872:38-160872:55',
	'160873:1-160873:147',
	'160874:1-160874:51',
	'160874:97-160874:113',
	'160939:1-160939:123',
	'160940:1-160940:79',
	'160942:1-160942:12',
	'160943:1-160943:54',
	'160955:1-160955:130',
	'160955:133-160955:138',
	'160955:140-160955:151',
	'160955:153-160955:154',
	'160955:156-160955:172',
	'160955:175-160955:201',
	'160955:204-160955:206',
	'160956:2-160956:65',
	'160957:1-160957:953',
	'160998:2-160998:252',
	'161008:2-161008:77',
	'161016:2-161016:300',
	'162762:1-162762:102',
	'162765:1-162765:40',
	'162803:60-162803:124',
	'162803:135-162803:139',
	'162808:1-162808:51',
	'162811:1-162811:340',
	'162822:73-162822:307',
	'162825:1-162825:184',
	'162826:1-162826:24',
	'162828:1-162828:85',
	'162909:54-162909:290',
	'163046:1-163046:133',
	'163046:135-163046:238',
	'163069:73-163069:452',
	'163069:468-163069:633',
	'163071:1-163071:161',
	'163078:1-163078:23',
	'163232:110-163232:149',
	'163233:1-163233:283',
	'163234:1-163234:66',
	'163235:1-163235:461',
	'163237:1-163237:213',
	'163238:9-163238:15',
	'163252:60-163252:137',
	'163255:1-163255:359',
	'163255:412-163255:844',
	'163255:846-163255:846',
	'163255:848-163255:977',
	'163261:1-163261:3',
	'163261:10-163261:126',
	'163270:1-163270:76',
	'163270:79-163270:96',
	'163270:99-163270:475',
	'163270:479-163270:527',
	'163270:529-163270:685',
	'163270:695-163270:928',
	'163286:112-163286:401',
	'163289:1-163289:388',
	'163296:59-163296:230',
	'163296:232-163296:585',
	'163297:1-163297:219',
	'163300:1-163300:616',
	'163301:1-163301:192',
	'163302:1-163302:190',
	'163332:43-163332:118',
	'163332:224-163332:264',
	'163332:266-163332:599',
	'163332:601-163332:639',
	'163332:641-163332:801',
	'163333:1-163333:106',
	'163334:1-163334:35',
	'163334:37-163334:37',
	'163334:156-163334:556',
	'163337:1-163337:18',
	'163337:27-163337:201',
	'163337:203-163337:426',
	'163337:434-163337:461',
	'163338:1-163338:164',
	'163339:1-163339:172',
	'163340:1-163340:488',
	'163358:39-163358:63',
	'163369:1-163369:94',
	'163370:1-163370:147',
	'163371:1-163371:107',
	'163371:148-163371:356',
	'163372:1-163372:52',
	'163374:1-163374:599',
	'163374:603-163374:863',
	'163375:1-163375:10',
	'163376:1-163376:20',
	'163376:22-163376:246',
	'163378:1-163378:81',
	'163378:89-163378:179',
	'163378:181-163378:272',
	'163378:306-163378:615',
	'163385:52-163385:240',
	'163385:244-163385:405',
	'163387:1-163387:256',
	'163402:37-163402:250',
	'163402:271-163402:560',
	'163402:581-163402:582',
	'163402:586-163402:801',
	'163475:30-163475:295',
	'163476:1-163476:94',
	'163476:98-163476:212',
	'163478:1-163478:70',
	'163479:1-163479:175',
	'163480:1-163480:92',
	'163480:96-163480:188',
	'163480:190-163480:191',
	'163481:1-163481:72',
	'163481:74-163481:77',
	'163481:79-163481:79',
	'163482:1-163482:27',
	'163482:48-163482:48',
	'163483:1-163483:57',
	'163582:1-163582:22',
	'163583:1-163583:63',
	'163583:65-163583:92',
	'163583:96-163583:155',
	'163583:157-163583:173',
	'163583:175-163583:219',
	'163584:1-163584:56',
	'163585:1-163585:32',
	'163586:1-163586:75',
	'163587:1-163587:52',
	'163588:1-163588:8',
	'163588:10-163588:446',
	'163589:1-163589:49',
	'163589:51-163589:160',
	'163596:1-163596:29',
	'163630:76-163630:164',
	'163630:176-163630:185',
	'163655:15-163655:23',
	'163657:1-163657:140',
	'163658:1-163658:3',
	'163659:1-163659:374',
	'163659:376-163659:650',
	'163659:652-163659:705',
	'163660:1-163660:74',
	'163661:1-163661:17',
	'163662:1-163662:154',
	'163663:1-163663:106',
	'163663:109-163663:246',
	'163664:1-163664:119',
	'163664:121-163664:178',
	'163668:1-163668:53',
	'163668:57-163668:136',
	'163668:140-163668:213',
	'163738:34-163738:311',
	'163757:1-163757:40',
	'163758:1-163758:17',
	'163758:19-163758:220',
	'163758:222-163758:224',
	'163758:236-163758:276',
	'163758:283-163758:374',
	'163758:376-163758:466',
	'163758:468-163758:591',
	'163759:1-163759:60',
	'163759:62-163759:72',
	'163759:74-163759:456',
	'163759:458-163759:461',
	'163759:463-163759:482',
	'163759:504-163759:510',
	'163760:1-163760:162',
	'163760:165-163760:340',
	'163761:1-163761:203',
	'163763:1-163763:79',
	'163765:1-163765:321',
	'163795:10-163795:34',
	'163795:36-163795:36',
	'163795:38-163795:43',
	'163796:1-163796:182',
	'163817:50-163817:140',
	'163817:154-163817:205',
	'163817:216-163817:295',
	'163817:305-163817:346',
	'163817:358-163817:457',
	'163817:561-163817:603',
	'163817:618-163817:966',
	'163869:79-163869:123',
	'165088:107-165088:266',
	'165098:124-165098:187',
	'165098:190-165098:193',
	'165098:195-165098:248',
	'165098:250-165098:254',
	'165098:256-165098:331',
	'165098:333-165098:367',
	'165098:369-165098:415',
	'165099:1-165099:105',
	'165102:1-165102:185',
	'165103:1-165103:440',
	'165120:82-165120:97',
	'165121:1-165121:466',
	'165205:80-165205:248',
	'165208:1-165208:101',
	'165364:45-165364:111',
	'165364:114-165364:147',
	'165364:160-165364:807',
	'165364:809-165364:1220',
	'165364:1260-165364:1301',
	'165402:1-165402:28',
	'165415:58-165415:85',
	'165415:88-165415:640',
	'165415:643-165415:707',
	'165415:712-165415:777',
	'165415:780-165415:1356',
	'165415:1360-165415:1383',
	'165467:39-165467:708',
	'165472:1-165472:184',
	'165472:186-165472:882',
	'165486:37-165486:102',
	'165487:1-165487:151',
	'165506:54-165506:170',
	'165514:72-165514:244',
	'165514:246-165514:270',
	'165514:283-165514:560',
	'165514:562-165514:567',
	'165548:1-165548:363',
	'165548:365-165548:381',
	'165548:384-165548:589',
	'165558:1-165558:62',
	'165567:54-165567:109',
	'165567:114-165567:309',
	'165567:315-165567:631',
	'165570:1-165570:2',
	'165570:5-165570:83',
	'165570:88-165570:207',
	'165570:209-165570:351',
	'165570:355-165570:942',
	'165570:944-165570:946',
	'165617:26-165617:52',
	'165617:54-165617:143',
	'165617:145-165617:288',
	'165620:14-165620:19',
	'165633:56-165633:62',
	'165633:64-165633:64',
	'165633:66-165633:129',
	'165633:133-165633:317',
	'165633:319-165633:500',
	'165970:67-165970:263',
	'165970:266-165970:329',
	'165970:331-165970:335',
	'165993:71-165993:873',
	'165993:879-165993:1660',
	'165993:1665-165993:1697',
	'166011:1-166011:81',
	'166011:83-166011:83',
	'166033:35-166033:53',
	'166033:59-166033:330',
	'166033:336-166033:355',
	'166033:360-166033:444',
	'166033:450-166033:606',
	'166033:613-166033:707',
	'166033:713-166033:1233',
	'166034:1-166034:109',
	'166034:115-166034:228',
	'166034:234-166034:307',
	'166049:53-166049:86',
	'166049:88-166049:236',
	'166049:242-166049:674',
	'166149:1-166149:2',
	'166150:1-166150:99',
	'166150:101-166150:116',
	'166161:38-166161:120',
	'166161:122-166161:123',
	'166161:126-166161:126',
	'166163:1-166163:12',
	'166163:14-166163:33',
	'166164:1-166164:32',
	'166164:34-166164:40',
	'166346:48-166346:210',
	'166346:212-166346:215',
	'166374:46-166374:64',
	'166374:66-166374:188',
	'166380:1-166380:367',
	'166380:373-166380:711',
	'166380:715-166380:1400',
	'166380:1406-166380:1809',
	'166408:67-166408:283',
	'166408:291-166408:947',
	'166408:953-166408:1236',
	'166429:33-166429:89',
	'166438:32-166438:85',
	'166438:87-166438:856',
	'166438:858-166438:866',
	'166462:78-166462:102',
	'166462:108-166462:317',
	'166462:323-166462:526',
	'166486:54-166486:75',
	'166486:80-166486:95',
	'166486:97-166486:174',
	'166502:43-166502:78',
	'166502:83-166502:109',
	'166512:42-166512:430',
	'166512:432-166512:487',
	'166512:491-166512:605',
	'166512:611-166512:1279',
	'166512:1281-166512:1818',
	'166512:1821-166512:1862',
	'166512:1868-166512:1868',
	'166512:1870-166512:1871',
	'166512:1873-166512:1874',
	'166514:1-166514:455',
	'166514:460-166514:464',
	'166530:43-166530:63',
	'166554:46-166554:218',
	'166554:224-166554:287',
	'166554:290-166554:317',
	'166554:320-166554:595',
	'166554:597-166554:730',
	'166554:732-166554:734',
	'166554:736-166554:736',
	'166563:1-166563:276',
	'166563:492-166563:748',
	'166565:1-166565:147',
	'166565:153-166565:312',
	'166565:316-166565:467',
	'166565:469-166565:898',
	'166699:55-166699:234',
	'166699:240-166699:415',
	'166699:421-166699:477',
	'166699:483-166699:677',
	'166699:681-166699:912',
	'166701:1-166701:13',
	'166701:16-166701:319',
	'166701:324-166701:506',
	'166701:513-166701:551',
	'166701:557-166701:672',
	'166701:681-166701:705',
	'166701:712-166701:724',
	'166701:731-166701:757',
	'166701:764-166701:777',
	'166701:783-166701:792',
	'166763:46-166763:168',
	'166763:174-166763:650',
	'166781:41-166781:111',
	'166781:115-166781:115',
	'166781:117-166781:233',
	'166781:236-166781:253',
	'166781:255-166781:382',
	'166782:1-166782:569',
	'166784:1-166784:114',
	'166784:119-166784:276',
	'166784:281-166784:365',
	'166787:1-166787:55',
	'166787:60-166787:127',
	'166787:132-166787:364',
	'166839:43-166839:173',
	'166839:178-166839:297',
	'166839:299-166839:302',
	'166841:1-166841:845',
	'166841:851-166841:876',
	'166841:882-166841:977',
	'166841:984-166841:984',
	'166841:988-166841:992',
	'166841:998-166841:1015',
	'166842:1-166842:170',
	'166859:62-166859:418',
	'166859:421-166859:421',
	'166859:423-166859:423',
	'166860:1-166860:21',
	'166861:1-166861:6',
	'166861:8-166861:13',
	'166864:1-166864:29',
	'166864:31-166864:77',
	'166864:79-166864:99',
	'166864:102-166864:119',
	'166864:125-166864:247',
	'166864:249-166864:307',
	'166864:311-166864:365',
	'166864:367-166864:374',
	'166864:378-166864:454',
	'166864:478-166864:537',
	'166888:56-166888:90',
	'166888:93-166888:154',
	'166888:156-166888:394',
	'166888:398-166888:470',
	'166889:1-166889:73',
	'166889:79-166889:228',
	'166890:1-166890:441',
	'166894:1-166894:190',
	'166895:1-166895:66',
	'166895:72-166895:597',
	'166895:599-166895:603',
	'166911:58-166911:76',
	'166911:81-166911:104',
	'166922:1-166922:39',
	'166922:41-166922:105',
	'166922:110-166922:340',
	'166922:345-166922:418',
	'166922:423-166922:747',
	'166922:752-166922:769',
	'166922:773-166922:773',
	'166923:1-166923:382',
	'166923:389-166923:470',
	'166946:41-166946:72',
	'166946:75-166946:201',
	'166950:1-166950:1',
	'166950:8-166950:31',
	'166950:36-166950:210',
	'166950:216-166950:877',
	'166950:883-166950:950',
	'166950:956-166950:1012',
	'166950:1018-166950:1321',
	'166950:1327-166950:1345',
	'166950:1347-166950:1438',
	'166960:1-166960:137',
	'166960:143-166960:166',
	'166966:1-166966:238',
	'166967:1-166967:220',
	'167039:20-167039:92',
	'167039:98-167039:228',
	'167041:1-167041:336',
	'167041:339-167041:391',
	'167041:396-167041:462',
	'167041:467-167041:663',
	'167043:1-167043:125',
	'167043:130-167043:235',
	'167078:40-167078:174',
	'167098:62-167098:90',
	'167098:92-167098:162',
	'167098:167-167098:406',
	'167098:448-167098:461',
	'167102:1-167102:42',
	'167102:48-167102:233',
	'167102:235-167102:317',
	'167102:323-167102:430',
	'167103:1-167103:94',
	'167151:1-167151:42',
	'167281:18-167281:140',
	'167281:146-167281:315',
	'167281:317-167281:593',
	'167282:1-167282:441',
	'167284:1-167284:315',
	'167284:320-167284:346',
	'167284:356-167284:395',
	'167284:399-167284:474',
	'167284:476-167284:1157',
	'167284:1160-167284:1628',
	'167284:1633-167284:1644',
	'167551:56-167551:190',
	'167551:196-167551:471',
	'167673:210-167673:236',
	'167673:239-167673:305',
	'167673:309-167673:418',
	'167673:423-167673:447',
	'167674:1-167674:345',
	'167675:1-167675:129',
	'167675:133-167675:299',
	'167675:301-167675:617',
	'167675:690-167675:707',
	'167675:710-167675:712',
	'167675:715-167675:716',
	'167675:719-167675:719',
	'167675:721-167675:725',
	'167675:740-167675:741',
	'167675:748-167675:758',
	'167675:762-167675:770',
	'167675:774-167675:787',
	'167675:793-167675:797',
	'167675:811-167675:1062',
	'167676:1-167676:278',
	'167676:289-167676:450',
	'167740:79-167740:126',
	'167740:132-167740:168',
	'167740:170-167740:173',
	'167746:56-167746:384',
	'167754:62-167754:103',
	'167784:51-167784:67',
	'167786:1-167786:1',
	'167786:11-167786:75',
	'167786:81-167786:176',
	'167807:60-167807:159',
	'167807:178-167807:204',
	'167807:210-167807:482',
	'167807:484-167807:558',
	'167807:560-167807:872',
	'167807:878-167807:1441',
	'167807:1444-167807:1842',
	'167830:1-167830:437',
	'167830:442-167830:587',
	'167830:590-167830:828',
	'167830:834-167830:1242',
	'167898:108-167898:619',
	'167898:621-167898:995',
	'167898:1001-167898:1010',
	'167898:1013-167898:1053',
	'167898:1057-167898:1295',
	'167898:1298-167898:1762',
	'167913:1-167913:126',
	'167913:128-167913:432',
        #--- August 5 re-reco ---#
	'170826:50-170826:122',
	'170826:139-170826:243',
	'170826:248-170826:308',
	'170842:1-170842:27',
	'170842:32-170842:96',
	'170842:102-170842:331',
	'170854:1-170854:336',
	'170854:341-170854:414',
	'170854:420-170854:470',
	'170854:475-170854:578',
	'170876:1-170876:110',
	'170876:115-170876:295',
	'170876:301-170876:516',
	'170876:518-170876:518',
	'170896:1-170896:212',
	'170899:1-170899:84',
	'170901:1-170901:153',
	'170901:159-170901:200',
	'171050:47-171050:74',
	'171050:80-171050:337',
	'171050:342-171050:369',
	'171050:371-171050:379',
	'171050:384-171050:423',
	'171050:427-171050:467',
	'171050:471-171050:648',
	'171091:1-171091:135',
	'171098:1-171098:8',
	'171102:1-171102:19',
	'171106:1-171106:27',
	'171106:32-171106:288',
	'171117:1-171117:54',
	'171117:56-171117:78',
	'171117:80-171117:84',
	'171156:42-171156:106',
	'171156:111-171156:686',
	'171156:688-171156:692',
	'171178:1-171178:92',
	'171178:97-171178:205',
	'171178:210-171178:551',
	'171178:556-171178:557',
	'171178:563-171178:787',
	'171178:792-171178:1043',
	'171219:48-171219:151',
	'171219:153-171219:162',
	'171274:88-171274:137',
	'171274:140-171274:143',
	'171282:1-171282:12',
	'171282:14-171282:99',
	'171282:104-171282:134',
	'171282:140-171282:171',
	'171315:53-171315:225',
	'171369:42-171369:130',
	'171369:136-171369:137',
	'171369:144-171369:161',
	'171446:58-171446:144',
	'171484:61-171484:202',
	'171484:207-171484:370',
	'171484:377-171484:432',
	'171578:47-171578:150',
	'171578:156-171578:174',
	'171578:179-171578:314',
	'171578:320-171578:347',
	'171578:353-171578:480',
	'171578:487-171578:572',
	'171578:578-171578:816',
	'171578:821-171578:974',
	'171812:59-171812:296',
	'171812:301-171812:438',
	'171876:1-171876:362',
	'171876:368-171876:391',
	'171876:397-171876:533',
	'171880:19-171880:202',
	'171895:34-171895:34',
	'171897:1-171897:23',
	'171897:205-171897:437',
	'171897:442-171897:511',
	'171897:517-171897:542',
	'171921:51-171921:141',
	'171926:1-171926:49',
	'171926:51-171926:155',
	'171926:161-171926:172',
	'171926:177-171926:264',
	'172014:1-172014:64',
	'172014:66-172014:143',
	'172014:149-172014:243',
	'172024:1-172024:46',
	'172033:1-172033:65',
	'172033:71-172033:277',
	'172033:282-172033:375',
	'172033:380-172033:473',
	'172033:478-172033:749',
	'172163:36-172163:109',
	'172163:115-172163:469',
	'172208:61-172208:199',
	'172252:32-172252:54',
	'172254:1-172254:35',
	'172255:1-172255:40',
	'172268:56-172268:169',
	'172286:52-172286:177',
	'172286:184-172286:216',
	'172389:34-172389:144',
	'172389:150-172389:428',
	'172389:434-172389:460',
	'172399:57-172399:226',
	'172400:1-172400:495',
	'172400:500-172400:690',
	'172400:696-172400:704',
	'172401:1-172401:2',
	'172401:5-172401:148',
	'172411:85-172411:349',
	'172478:1-172478:110',
	'172619:1-172619:77',
        #--- Prompt V6 ---#
	'172620:1-172620:495',
	'172630:36-172630:64',
	'172630:68-172630:134',
	'172630:139-172630:160',
	'172635:1-172635:18',
	'172635:24-172635:268',
	'172778:48-172778:97',
	'172791:65-172791:413',
	'172791:418-172791:569',
	'172791:571-172791:715',
	'172791:721-172791:1295',
	'172791:1300-172791:1536',
	'172791:1542-172791:1645',
	'172791:1649-172791:1658',
	'172798:1-172798:31',
	'172799:1-172799:367',
	'172801:1-172801:679',
	'172801:681-172801:750',
	'172801:753-172801:766',
	'172801:768-172801:815',
	'172801:819-172801:837',
	'172801:839-172801:861',
	'172801:863-172801:909',
	'172801:911-172801:1139',
	'172802:1-172802:629',
	'172802:634-172802:784',
	'172819:57-172819:87',
	'172819:92-172819:254',
	'172822:1-172822:75',
	'172822:77-172822:127',
	'172822:133-172822:353',
	'172822:358-172822:662',
	'172822:667-172822:832',
	'172822:837-172822:1096',
	'172822:1102-172822:1112',
	'172822:1118-172822:2121',
	'172822:2127-172822:2333',
	'172824:1-172824:54',
	'172847:62-172847:131',
	'172865:36-172865:151',
	'172865:157-172865:382',
	'172868:1-172868:29',
	'172868:34-172868:708',
	'172868:714-172868:909',
	'172868:912-172868:1829',
	'172868:1835-172868:1970',
	'172949:55-172949:127',
	'172949:133-172949:159',
	'172949:161-172949:263',
	'172949:269-172949:274',
	'172949:276-172949:892',
	'172949:898-172949:928',
	'172949:933-172949:1287',
	'172951:1-172951:52',
	'172952:1-172952:195',
	'172952:197-172952:670',
	'172952:675-172952:797',
	'172952:799-172952:1164',
	'172952:1169-172952:1562',
	'172992:505-172992:511',
	'172992:516-172992:713',
	'172992:752-172992:874',
	'172992:880-172992:946',
	'172999:1-172999:71',
	'172999:77-172999:294',
	'173198:49-173198:322',
	'173198:324-173198:356',
	'173198:362-173198:485',
	'173198:488-173198:768',
	'173198:770-173198:823'
])
