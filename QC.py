import pandas as pd

bowtie_df = pd.read_csv('./GEO/bowtie_srr_results.csv')
bowtie_df = bowtie_df.set_index(['sample_id'])
NSC_RSC_df = pd.read_csv('./GEO/NSC_RSC_SRR_sample2.txt', sep='\t')
NSC_RSC_df['Run_ID'] = NSC_RSC_df['Run_ID'].str.replace('.tagAlign', '')
NSC_RSC_df = NSC_RSC_df.set_index(['Run_ID'])

NRF_df = pd.read_csv('./GEO/NRF_GEO.csv')
NRF_df = NRF_df.set_index(['sample_name'])

pair_df = pd.read_csv('./GEO/GEO_danpos_pair.csv')
pair_df = pair_df.set_index(['Run_ID'])

final_df = pair_df.join(bowtie_df).join(NSC_RSC_df).join(NRF_df)

final_df.to_csv('GEO_QC_final.csv')

a = ['GSM1003709', 'GSM1003712', 'GSM1003715', 'GSM1003718', 'GSM1003721', 'GSM1003724', 'GSM1003729', 'GSM1003732', 'GSM1003735', 'GSM1003738', 'GSM1003741', 'GSM1006722', 'GSM1007996', 'GSM1020362', 'GSM1033552', 'GSM1033553', 'GSM1033554', 'GSM1033555', 'GSM1033556', 'GSM1033557', 'GSM1033570', 'GSM1033571', 'GSM1033572', 'GSM1033573', 'GSM1033574', 'GSM1033575', 'GSM1033758', 'GSM1035429', 'GSM1035438', 'GSM1055816', 'GSM1055817', 'GSM1065009', 'GSM1067280', 'GSM1080936', 'GSM1089821', 'GSM1122656', 'GSM1122657', 'GSM1122666', 'GSM1129348', 'GSM1129350', 'GSM1129352', 'GSM1129354', 'GSM1133648', 'GSM1133649', 'GSM1135046', 'GSM1135047', 'GSM1141672', 'GSM1146439', 'GSM1146440', 'GSM1146441', 'GSM1146442', 'GSM1146450', 'GSM1203228', 'GSM1204472', 'GSM1204473', 'GSM1208811', 'GSM1217919', 'GSM1217922', 'GSM1217924', 'GSM1217925', 'GSM1217926', 'GSM1217928', 'GSM1217929', 'GSM1217930', 'GSM1217931', 'GSM1217952', 'GSM1233880', 'GSM1233881', 'GSM1233905', 'GSM1233906', 'GSM1233907', 'GSM1233926', 'GSM1233927', 'GSM1233947', 'GSM1233948', 'GSM1233949', 'GSM1233969', 'GSM1233970', 'GSM1233971', 'GSM1233988', 'GSM1233989', 'GSM1234003', 'GSM1234004', 'GSM1234020', 'GSM1234021', 'GSM1234037', 'GSM1234038', 'GSM1234054', 'GSM1234055', 'GSM1234071', 'GSM1234072', 'GSM1234091', 'GSM1234092', 'GSM1234093', 'GSM1234113', 'GSM1234114', 'GSM1234115', 'GSM1234135', 'GSM1234136', 'GSM1234137', 'GSM1234138', 'GSM1234155', 'GSM1234156', 'GSM1234173', 'GSM1234174', 'GSM1234191', 'GSM1234192', 'GSM1234209', 'GSM1234210', 'GSM1234232', 'GSM1234233', 'GSM1240115', 'GSM1240116', 'GSM1261675', 'GSM1262399', 'GSM1262400', 'GSM1265857', 'GSM1265858', 'GSM1273644', 'GSM1273649', 'GSM1273654', 'GSM1273659', 'GSM1277990', 'GSM1295085', 'GSM1295095', 'GSM1295099', 'GSM1303758', 'GSM1305207', 'GSM1305208', 'GSM1305209', 'GSM1313524', 'GSM1313527', 'GSM1320313', 'GSM1320317', 'GSM1320321', 'GSM1320325', 'GSM1322273', 'GSM1325254', 'GSM1325255', 'GSM1326437', 'GSM1326438', 'GSM1326439', 'GSM1326440', 'GSM1326441', 'GSM1326442', 'GSM1338273', 'GSM1338274', 'GSM1338275', 'GSM1338285', 'GSM1338286', 'GSM1338287', 'GSM1354430', 'GSM1354431', 'GSM1356566', 'GSM1359503', 'GSM1365901', 'GSM1372853', 'GSM1372859', 'GSM1372864', 'GSM1372871', 'GSM1382470', 'GSM1382471', 'GSM1382480', 'GSM1382481', 'GSM1383856', 'GSM1383862', 'GSM1383868', 'GSM1383874', 'GSM1383881', 'GSM1385784', 'GSM1385787', 'GSM1386345', 'GSM1402463', 'GSM1410781', 'GSM1415874', 'GSM1415882', 'GSM1418959', 'GSM1418960', 'GSM1420159', 'GSM1436881', 'GSM1436882', 'GSM1447335', 'GSM1462465', 'GSM1462466', 'GSM1479215', 'GSM1482819', 'GSM1482821', 'GSM1482878', 'GSM1486000', 'GSM1486001', 'GSM1498896', 'GSM1498897', 'GSM1498898', 'GSM1498899', 'GSM1499823', 'GSM1499824', 'GSM1501176', 'GSM1508939', 'GSM1508945', 'GSM1513832', 'GSM1513833', 'GSM1513893', 'GSM1513894', 'GSM1513897', 'GSM1513898', 'GSM1519153', 'GSM1519154', 'GSM1519155', 'GSM1519156', 'GSM1519157', 'GSM1519158', 'GSM1519159', 'GSM1519160', 'GSM1519161', 'GSM1519162', 'GSM1519163', 'GSM1519164', 'GSM1519165', 'GSM1519166', 'GSM1519167', 'GSM1519168', 'GSM1519169', 'GSM1519170', 'GSM1519171', 'GSM1519172', 'GSM1519173', 'GSM1519174', 'GSM1519175', 'GSM1519176', 'GSM1519177', 'GSM1519178', 'GSM1519179', 'GSM1519180', 'GSM1519181', 'GSM1519182', 'GSM1519183', 'GSM1519184', 'GSM1519185', 'GSM1519186', 'GSM1519187', 'GSM1519188', 'GSM1519189', 'GSM1519190', 'GSM1519191', 'GSM1519192', 'GSM1519193', 'GSM1519194', 'GSM1519195', 'GSM1519196', 'GSM1519197', 'GSM1519198', 'GSM1521724', 'GSM1521729', 'GSM1521733', 'GSM1521738', 'GSM1521743', 'GSM1521748', 'GSM1521753', 'GSM1521758', 'GSM1527677', 'GSM1527679', 'GSM1527826', 'GSM1527827', 'GSM1534445', 'GSM1536034', 'GSM1536038', 'GSM1537283', 'GSM1537284', 'GSM1537285', 'GSM1537286', 'GSM1537287', 'GSM1541017', 'GSM1541018', 'GSM1552410', 'GSM1552413', 'GSM1556824', 'GSM1556825', 'GSM1565791', 'GSM1565798', 'GSM1565805', 'GSM1565812', 'GSM1568245', 'GSM1569565', 'GSM1569573', 'GSM1569577', 'GSM1572767', 'GSM1572770', 'GSM1572771', 'GSM1572774', 'GSM1572775', 'GSM1572777', 'GSM1572778', 'GSM1572779', 'GSM1572780', 'GSM1572781', 'GSM1572783', 'GSM1572784', 'GSM1572785', 'GSM1572786', 'GSM1572787', 'GSM1572788', 'GSM1572789', 'GSM1574256', 'GSM1574261', 'GSM1599147', 'GSM1599148', 'GSM1599149', 'GSM1599150', 'GSM1599151', 'GSM1599152', 'GSM1599153', 'GSM1599154', 'GSM1602241', 'GSM1607388', 'GSM1607397', 'GSM1607401', 'GSM1607402', 'GSM1607403', 'GSM1613314', 'GSM1615880', 'GSM1615881', 'GSM1647368', 'GSM1647369', 'GSM1647370', 'GSM1647371', 'GSM1648034', 'GSM1648037', 'GSM1653247', 'GSM1653249', 'GSM1653251', 'GSM1663091', 'GSM1663092', 'GSM1663094', 'GSM1663096', 'GSM1666202', 'GSM1666203', 'GSM1666204', 'GSM1666205', 'GSM1666384', 'GSM1668440', 'GSM1668443', 'GSM1668446', 'GSM1670701', 'GSM1688583', 'GSM1688584', 'GSM1688585', 'GSM1688586', 'GSM1688587', 'GSM1688588', 'GSM1688590', 'GSM1688591', 'GSM1688592', 'GSM1688593', 'GSM1688594', 'GSM1688595', 'GSM1697670', 'GSM1697673', 'GSM1697676', 'GSM1700380', 'GSM1700381', 'GSM1700386', 'GSM1700387', 'GSM1700392', 'GSM1700393', 'GSM1708663', 'GSM1717082', 'GSM1717083', 'GSM1722606', 'GSM1722607', 'GSM1722608', 'GSM1722609', 'GSM1782766', 'GSM1816070', 'GSM1816072', 'GSM1816684', 'GSM1816685', 'GSM1817174', 'GSM1817175', 'GSM1817176', 'GSM1821389', 'GSM1824910', 'GSM1824911', 'GSM1832644', 'GSM1832645', 'GSM1835837', 'GSM1835840', 'GSM1835843', 'GSM1835846', 'GSM1835849', 'GSM1835852', 'GSM1835859', 'GSM1835875', 'GSM1835880', 'GSM1835886', 'GSM1835889', 'GSM1835894', 'GSM1835899', 'GSM1835986', 'GSM1835988', 'GSM1841291', 'GSM1841302', 'GSM1846779', 'GSM1846780', 'GSM1846782', 'GSM1846784', 'GSM1846800', 'GSM1846804', 'GSM1847924', 'GSM1847925', 'GSM1847926', 'GSM1847927', 'GSM1847928', 'GSM1847929', 'GSM1847930', 'GSM1847931', 'GSM1847932', 'GSM1847933', 'GSM1847934', 'GSM1847935', 'GSM1847936', 'GSM1847937', 'GSM1847938', 'GSM1847939', 'GSM1847940', 'GSM1847941', 'GSM1847942', 'GSM1847943', 'GSM1847944', 'GSM1847945', 'GSM1847946', 'GSM1847947', 'GSM1847948', 'GSM1847949', 'GSM1847950', 'GSM1847951', 'GSM1847952', 'GSM1847953', 'GSM1847954', 'GSM1847955', 'GSM1847956', 'GSM1847957', 'GSM1847958', 'GSM1847959', 'GSM1847960', 'GSM1847961', 'GSM1847962', 'GSM1847963', 'GSM1847964', 'GSM1847965', 'GSM1847966', 'GSM1847967', 'GSM1847968', 'GSM1847969', 'GSM1847970', 'GSM1847971', 'GSM1847972', 'GSM1847973', 'GSM1847974', 'GSM1847975', 'GSM1847976', 'GSM1847977', 'GSM1847978', 'GSM1847979', 'GSM1847980', 'GSM1847981', 'GSM1847982', 'GSM1847983', 'GSM1853822', 'GSM1853824', 'GSM1855843', 'GSM1855844', 'GSM1855845', 'GSM1855846', 'GSM1863071', 'GSM1863074', 'GSM1863077', 'GSM1866053', 'GSM1866056', 'GSM1866086', 'GSM1866092', 'GSM1866096', 'GSM1866100', 'GSM1866117', 'GSM1866119', 'GSM1866120', 'GSM1866121', 'GSM1866125', 'GSM1866127', 'GSM1866129', 'GSM1866131', 'GSM1866133', 'GSM1866135', 'GSM1866137', 'GSM1866699', 'GSM1868870', 'GSM1868871', 'GSM1869134', 'GSM1869140', 'GSM1869146', 'GSM1871944', 'GSM1871945', 'GSM1871946', 'GSM1871947', 'GSM1871948', 'GSM1871949', 'GSM1874929', 'GSM1874930', 'GSM1874931', 'GSM1888683', 'GSM1888685', 'GSM1888687', 'GSM1888689', 'GSM1888691', 'GSM1888693', 'GSM1888695', 'GSM1888697', 'GSM1888699', 'GSM1888701', 'GSM1888703', 'GSM1888705', 'GSM1888707', 'GSM1888709', 'GSM1888711', 'GSM1888713', 'GSM1888732', 'GSM1888736', 'GSM1888740', 'GSM1888744', 'GSM1888748', 'GSM1888752', 'GSM1888756', 'GSM1888760', 'GSM1888763', 'GSM1888765', 'GSM1888767', 'GSM1888769', 'GSM1888771', 'GSM1888773', 'GSM1888775', 'GSM1888777', 'GSM1889909', 'GSM1889910', 'GSM1892605', 'GSM1892653', 'GSM1918513', 'GSM1918514', 'GSM1918515', 'GSM1918516', 'GSM1918517', 'GSM1918518', 'GSM1918519', 'GSM1918520', 'GSM1918521', 'GSM1918522', 'GSM1918523', 'GSM1918524', 'GSM1918525', 'GSM1918526', 'GSM1918527', 'GSM1918528', 'GSM1918529', 'GSM1918530', 'GSM1918531', 'GSM1918532', 'GSM1918533', 'GSM1918534', 'GSM1918535', 'GSM1918536', 'GSM1918597', 'GSM1918598', 'GSM1918599', 'GSM1918600', 'GSM1918601', 'GSM1918602', 'GSM1918603', 'GSM1918604', 'GSM1918605', 'GSM1918606', 'GSM1918607', 'GSM1918608', 'GSM1918609', 'GSM1918610', 'GSM1918611', 'GSM1918612', 'GSM1918613', 'GSM1918614', 'GSM1918615', 'GSM1918616', 'GSM1918637', 'GSM1918638', 'GSM1918639', 'GSM1918640', 'GSM1918641', 'GSM1918642', 'GSM1918643', 'GSM1918644', 'GSM1918645', 'GSM1918646', 'GSM1918647', 'GSM1918648', 'GSM1918649', 'GSM1918650', 'GSM1918651', 'GSM1918652', 'GSM1918653', 'GSM1918654', 'GSM1918655', 'GSM1918656', 'GSM1922952', 'GSM1922953', 'GSM1922954', 'GSM1922955', 'GSM1922956', 'GSM1922957', 'GSM1922958', 'GSM1922959', 'GSM1922960', 'GSM1922961', 'GSM1934090', 'GSM1958042', 'GSM1958043', 'GSM1960245', 'GSM1960246', 'GSM1960247', 'GSM1960248', 'GSM1961542', 'GSM1976295', 'GSM1976305', 'GSM2029137', 'GSM2029138', 'GSM2029139', 'GSM2029140', 'GSM2029141', 'GSM2029142', 'GSM2029143', 'GSM2029144', 'GSM2029145', 'GSM2029146', 'GSM2029147', 'GSM2029350', 'GSM2029355', 'GSM2029360', 'GSM2029365', 'GSM2029585', 'GSM2029590', 'GSM2046861', 'GSM2046862', 'GSM2046863', 'GSM2046864', 'GSM2046865', 'GSM2046866', 'GSM2047024', 'GSM2048288', 'GSM2048289', 'GSM2048290', 'GSM2048295', 'GSM2048296', 'GSM2048297', 'GSM2048302', 'GSM2048303', 'GSM2048307', 'GSM2048308', 'GSM2048312', 'GSM2048313', 'GSM2058901', 'GSM2058902', 'GSM2058907', 'GSM2058908', 'GSM2058913', 'GSM2058914', 'GSM2060761', 'GSM2060762', 'GSM2066617', 'GSM2066618', 'GSM2066619', 'GSM2067930', 'GSM2067931', 'GSM2067932', 'GSM2068358', 'GSM2068366', 'GSM2072641', 'GSM2120699', 'GSM2120704', 'GSM2120711', 'GSM2127460', 'GSM2127466', 'GSM2131286', 'GSM2131287', 'GSM2140614', 'GSM2140615', 'GSM2140616', 'GSM2140617', 'GSM2140618', 'GSM2159749', 'GSM2159751', 'GSM2159753', 'GSM2166074', 'GSM2166075', 'GSM2171845', 'GSM2171846', 'GSM2171847', 'GSM2171848', 'GSM2171849', 'GSM2171850', 'GSM2171851', 'GSM2171852', 'GSM2175563', 'GSM2175564', 'GSM2175579', 'GSM2175580', 'GSM2178507', 'GSM2178509', 'GSM2178511', 'GSM2178513', 'GSM2178515', 'GSM2178517', 'GSM2178519', 'GSM2178521', 'GSM2178523', 'GSM2178525', 'GSM2178527', 'GSM2178529', 'GSM2187240', 'GSM2187253', 'GSM2187255', 'GSM2187257', 'GSM2187260', 'GSM2187262', 'GSM2187264', 'GSM2279960', 'GSM2279968', 'GSM2279978', 'GSM2279987', 'GSM2279997', 'GSM2280008', 'GSM2280016', 'GSM2280024', 'GSM2280031', 'GSM2280039', 'GSM2330571', 'GSM2330572', 'GSM2337943', 'GSM2337944', 'GSM2337951', 'GSM2337952', 'GSM2341637', 'GSM2341644', 'GSM2356337', 'GSM2356343', 'GSM2357601', 'GSM2357609', 'GSM2357617', 'GSM2410096', 'GSM2410097', 'GSM2410098', 'GSM2410099', 'GSM2422846', 'GSM2422847', 'GSM2422850', 'GSM2422851', 'GSM2442785', 'GSM569085', 'GSM586869', 'GSM586872', 'GSM586873', 'GSM586880', 'GSM586885', 'GSM593365', 'GSM593366', 'GSM593367', 'GSM602296', 'GSM602845', 'GSM641814', 'GSM641817', 'GSM644990', 'GSM648495', 'GSM686935', 'GSM714809', 'GSM721134', 'GSM727572', 'GSM727596', 'GSM732911', 'GSM773476', 'GSM779435', 'GSM790215', 'GSM865280', 'GSM865281', 'GSM865282', 'GSM883688', 'GSM883689', 'GSM883690', 'GSM883691', 'GSM883692', 'GSM883693', 'GSM894068', 'GSM894069', 'GSM894070', 'GSM894071', 'GSM894072', 'GSM894073', 'GSM894084', 'GSM916966', 'GSM916968', 'GSM916971', 'GSM916975', 'GSM916977', 'GSM916980', 'GSM922787', 'GSM941545', 'GSM942110', 'GSM942111', 'GSM942112', 'GSM942113', 'GSM942114', 'GSM942115', 'GSM942116', 'GSM942117', 'GSM942118', 'GSM950844', 'GSM950847', 'GSM950850', 'GSM952458', 'GSM952466', 'GSM955983', 'GSM955984', 'GSM969569', 'GSM984397', 'GSM984398', 'GSM984399', 'GSM984400', 'GSM995942']

# not_include = []
# for s in final_df['Data_ID']:
#     # print s
#     if s not in a:
#         print s, 'not in the wig'
#         not_include.append(s)
#
# for aa in a:
#     if aa not in final_df['Data_ID'].values:
#         print aa, 'wig not in the list'
#
# print not_include
final_df['Run_ID'] = final_df.index
final_df = final_df.set_index(['Data_ID'])
final_df = final_df.ix[a, :]

final_df.to_csv('GEO_QC_final.csv')