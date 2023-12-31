OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63859343) q[0];
sx q[0];
rz(-0.25829691) q[0];
sx q[0];
rz(-2.5845085) q[0];
rz(-pi) q[1];
rz(0.090473526) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(-2.5475516) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0253804) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(1.8616574) q[1];
rz(-pi) q[2];
x q[2];
rz(2.071322) q[3];
sx q[3];
rz(-2.4300017) q[3];
sx q[3];
rz(-0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(-0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(3.0965565) q[0];
rz(-2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(2.3388458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40753663) q[0];
sx q[0];
rz(-1.7965839) q[0];
sx q[0];
rz(-1.9682103) q[0];
rz(1.0674092) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(-0.78546055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3992577) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(-1.6595112) q[1];
rz(0.12985274) q[3];
sx q[3];
rz(-2.2359071) q[3];
sx q[3];
rz(-3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7654483) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1012672) q[0];
sx q[0];
rz(-0.87653941) q[0];
sx q[0];
rz(0.56157748) q[0];
rz(0.88111934) q[2];
sx q[2];
rz(-1.8124049) q[2];
sx q[2];
rz(-1.2348262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9078319) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(0.16252653) q[1];
rz(-0.6891059) q[3];
sx q[3];
rz(-1.8854685) q[3];
sx q[3];
rz(0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-0.21425042) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4857805) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(-0.73894545) q[0];
x q[1];
rz(-1.8604943) q[2];
sx q[2];
rz(-2.2640267) q[2];
sx q[2];
rz(-1.9300269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2506634) q[1];
sx q[1];
rz(-2.2729985) q[1];
sx q[1];
rz(0.82171792) q[1];
rz(-pi) q[2];
rz(-1.7188844) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9966853) q[0];
sx q[0];
rz(-0.7888182) q[0];
sx q[0];
rz(-0.44992723) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(1.8051403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0817464) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(-2.5043143) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(-1.050625) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-3.040722) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8773008) q[0];
sx q[0];
rz(-1.722324) q[0];
sx q[0];
rz(-0.63440462) q[0];
x q[1];
rz(0.48730536) q[2];
sx q[2];
rz(-1.8308182) q[2];
sx q[2];
rz(0.11815182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
x q[2];
rz(-2.6359667) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(-1.0491919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(-2.4538453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6763517) q[0];
sx q[0];
rz(-2.7382593) q[0];
sx q[0];
rz(-2.1562735) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96104788) q[2];
sx q[2];
rz(-2.7597722) q[2];
sx q[2];
rz(1.0605937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75377611) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(0.17863518) q[1];
x q[2];
rz(2.1391159) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872902) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(1.646423) q[0];
rz(1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(-1.5471293) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5115576) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(-0.10314421) q[1];
rz(-pi) q[2];
rz(-0.49097455) q[3];
sx q[3];
rz(-0.51530616) q[3];
sx q[3];
rz(-1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(2.2081597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.371884) q[0];
sx q[0];
rz(-3.0196307) q[0];
sx q[0];
rz(2.6321649) q[0];
rz(1.1912212) q[2];
sx q[2];
rz(-1.2588725) q[2];
sx q[2];
rz(1.003007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48298207) q[1];
sx q[1];
rz(-2.254199) q[1];
sx q[1];
rz(-1.618209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4890355) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(2.488234) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-2.9577589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18098022) q[0];
sx q[0];
rz(-1.6913337) q[0];
sx q[0];
rz(1.4535849) q[0];
rz(-pi) q[1];
rz(-2.3887563) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(2.2019049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(-0.88840719) q[1];
x q[2];
rz(2.5582383) q[3];
sx q[3];
rz(-1.5577941) q[3];
sx q[3];
rz(2.0875723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(0.72485483) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(0.98060676) q[3];
sx q[3];
rz(-1.7173613) q[3];
sx q[3];
rz(2.2955256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
