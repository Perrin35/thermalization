OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6440808) q[0];
sx q[0];
rz(-1.3389791) q[0];
sx q[0];
rz(-0.31565491) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(-2.5375836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8020129) q[0];
sx q[0];
rz(-0.30671841) q[0];
sx q[0];
rz(1.40236) q[0];
rz(-pi) q[1];
rz(0.96678218) q[2];
sx q[2];
rz(-2.3748715) q[2];
sx q[2];
rz(-1.3441835) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4872514) q[1];
sx q[1];
rz(-0.62955925) q[1];
sx q[1];
rz(2.1653324) q[1];
rz(1.2617926) q[3];
sx q[3];
rz(-1.6741535) q[3];
sx q[3];
rz(2.7637533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.6319298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5556521) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(-1.9441388) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(-1.7053568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2510808) q[0];
sx q[0];
rz(-2.5174826) q[0];
sx q[0];
rz(1.9735543) q[0];
x q[1];
rz(-0.70296944) q[2];
sx q[2];
rz(-1.6466993) q[2];
sx q[2];
rz(0.10026201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20138559) q[1];
sx q[1];
rz(-1.6734945) q[1];
sx q[1];
rz(-2.9443355) q[1];
rz(1.7800434) q[3];
sx q[3];
rz(-2.0179141) q[3];
sx q[3];
rz(-2.1802878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88840914) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-0.02296981) q[2];
rz(2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70812923) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(2.1500812) q[0];
rz(-2.1082711) q[1];
sx q[1];
rz(-2.8585377) q[1];
sx q[1];
rz(1.2352157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9231222) q[0];
sx q[0];
rz(-1.3231131) q[0];
sx q[0];
rz(2.5367303) q[0];
x q[1];
rz(-0.60977817) q[2];
sx q[2];
rz(-1.9710566) q[2];
sx q[2];
rz(-0.61851172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4146862) q[1];
sx q[1];
rz(-0.79830805) q[1];
sx q[1];
rz(1.9195791) q[1];
x q[2];
rz(0.97930538) q[3];
sx q[3];
rz(-2.0404173) q[3];
sx q[3];
rz(-2.8036433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.027355) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(-0.99119622) q[2];
rz(-1.8441955) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(1.7245002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7436413) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(0.81746307) q[0];
rz(0.3262597) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(-0.78525966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2366339) q[0];
sx q[0];
rz(-1.3044622) q[0];
sx q[0];
rz(-1.8840709) q[0];
rz(-1.3738858) q[2];
sx q[2];
rz(-2.2047538) q[2];
sx q[2];
rz(-1.9206573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7689159) q[1];
sx q[1];
rz(-1.1296037) q[1];
sx q[1];
rz(-2.6720474) q[1];
rz(-pi) q[2];
rz(-1.4912075) q[3];
sx q[3];
rz(-1.2168988) q[3];
sx q[3];
rz(0.45512629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12104812) q[2];
sx q[2];
rz(-2.497017) q[2];
sx q[2];
rz(2.5784967) q[2];
rz(2.2480887) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(-1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9017482) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(1.595994) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-2.0183759) q[1];
sx q[1];
rz(-2.9603069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.793042) q[0];
sx q[0];
rz(-0.53213813) q[0];
sx q[0];
rz(1.1842968) q[0];
rz(-2.6659545) q[2];
sx q[2];
rz(-1.5283094) q[2];
sx q[2];
rz(-1.3459537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5511032) q[1];
sx q[1];
rz(-0.22399513) q[1];
sx q[1];
rz(1.444054) q[1];
rz(-pi) q[2];
rz(-0.44293176) q[3];
sx q[3];
rz(-1.7128403) q[3];
sx q[3];
rz(2.1403109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1229317) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(1.2459416) q[2];
rz(2.5490226) q[3];
sx q[3];
rz(-1.6962681) q[3];
sx q[3];
rz(0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6084006) q[0];
sx q[0];
rz(-0.99219457) q[0];
sx q[0];
rz(2.8327508) q[0];
rz(2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(-3.0016532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788743) q[0];
sx q[0];
rz(-2.1669183) q[0];
sx q[0];
rz(-2.8403502) q[0];
rz(-pi) q[1];
rz(-1.762319) q[2];
sx q[2];
rz(-0.46326783) q[2];
sx q[2];
rz(-2.6216216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9819239) q[1];
sx q[1];
rz(-0.86729151) q[1];
sx q[1];
rz(-2.7419006) q[1];
x q[2];
rz(-1.1902963) q[3];
sx q[3];
rz(-0.69529136) q[3];
sx q[3];
rz(-2.8075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(0.24197401) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.022843) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(1.2710849) q[0];
rz(-0.56308833) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(-1.44106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1313167) q[0];
sx q[0];
rz(-0.97945853) q[0];
sx q[0];
rz(0.32846668) q[0];
x q[1];
rz(-2.8071324) q[2];
sx q[2];
rz(-0.40142599) q[2];
sx q[2];
rz(-0.73094598) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59380925) q[1];
sx q[1];
rz(-1.9873706) q[1];
sx q[1];
rz(-1.814247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8981179) q[3];
sx q[3];
rz(-2.4015275) q[3];
sx q[3];
rz(-2.8748663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21961221) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(1.8249576) q[2];
rz(-2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(-1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.104326) q[0];
sx q[0];
rz(-0.21793652) q[0];
sx q[0];
rz(0.88687801) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(1.0221457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048269317) q[0];
sx q[0];
rz(-1.5962068) q[0];
sx q[0];
rz(-1.4396458) q[0];
rz(0.81424197) q[2];
sx q[2];
rz(-0.98746429) q[2];
sx q[2];
rz(-0.63177201) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8274668) q[1];
sx q[1];
rz(-2.3955401) q[1];
sx q[1];
rz(-1.5681727) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5660146) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(0.23250599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48978051) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(-2.1583648) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(-2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1061358) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5863093) q[0];
rz(-0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(0.78479016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0138088) q[0];
sx q[0];
rz(-2.0079566) q[0];
sx q[0];
rz(0.60370914) q[0];
rz(-pi) q[1];
x q[1];
rz(1.773514) q[2];
sx q[2];
rz(-1.7200816) q[2];
sx q[2];
rz(-1.210807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49336772) q[1];
sx q[1];
rz(-2.7052212) q[1];
sx q[1];
rz(3.0619951) q[1];
rz(-pi) q[2];
rz(1.9905213) q[3];
sx q[3];
rz(-1.0015206) q[3];
sx q[3];
rz(-1.9939533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67184225) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78275457) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(0.77761039) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(0.75540677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841254) q[0];
sx q[0];
rz(-1.3608772) q[0];
sx q[0];
rz(0.0046878417) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.586953) q[2];
sx q[2];
rz(2.0153869) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3302119) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(2.5274171) q[1];
rz(-pi) q[2];
rz(1.2437263) q[3];
sx q[3];
rz(-0.9216412) q[3];
sx q[3];
rz(-0.4285194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5767673) q[2];
sx q[2];
rz(-2.8456523) q[2];
sx q[2];
rz(2.9887065) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.2234756) q[1];
sx q[1];
rz(-1.9203095) q[1];
sx q[1];
rz(-0.65345678) q[1];
rz(2.5443947) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(0.85930227) q[3];
sx q[3];
rz(-0.51325428) q[3];
sx q[3];
rz(-0.2556066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
