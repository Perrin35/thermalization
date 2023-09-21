OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(1.7095837) q[0];
x q[1];
rz(-1.4957501) q[2];
sx q[2];
rz(-1.4805761) q[2];
sx q[2];
rz(-2.1580632) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0253804) q[1];
sx q[1];
rz(-2.930021) q[1];
sx q[1];
rz(1.2799353) q[1];
rz(-pi) q[2];
rz(-2.071322) q[3];
sx q[3];
rz(-2.4300017) q[3];
sx q[3];
rz(-2.4430032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(-1.9682103) q[0];
rz(-1.0674092) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(-2.3561321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-0.90157585) q[3];
sx q[3];
rz(-1.4687317) q[3];
sx q[3];
rz(-1.4249742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(-2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1012672) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-2.5800152) q[0];
rz(-pi) q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(-0.14112976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50946402) q[1];
sx q[1];
rz(-2.8948935) q[1];
sx q[1];
rz(-2.2798666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1717103) q[3];
sx q[3];
rz(-2.2200924) q[3];
sx q[3];
rz(0.89023512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65581215) q[0];
sx q[0];
rz(-0.27944767) q[0];
sx q[0];
rz(-2.4026472) q[0];
x q[1];
rz(-2.8104086) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(0.77510288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22073332) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-3.1130476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.9821092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75393049) q[0];
sx q[0];
rz(-1.8844863) q[0];
sx q[0];
rz(0.73648209) q[0];
rz(-pi) q[1];
rz(2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(-1.3364524) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24134484) q[1];
sx q[1];
rz(-2.4240139) q[1];
sx q[1];
rz(2.248583) q[1];
rz(-pi) q[2];
rz(-2.106835) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(1.1353726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3374702) q[0];
sx q[0];
rz(-0.94479783) q[0];
sx q[0];
rz(-1.7581598) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-3.0234408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(0.15921758) q[1];
x q[2];
rz(0.50562596) q[3];
sx q[3];
rz(-0.45590948) q[3];
sx q[3];
rz(-2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7589384) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9822385) q[0];
sx q[0];
rz(-1.9040477) q[0];
sx q[0];
rz(-2.9100145) q[0];
x q[1];
rz(0.22600941) q[2];
sx q[2];
rz(-1.2603716) q[2];
sx q[2];
rz(1.7058536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2076599) q[1];
sx q[1];
rz(-1.4414409) q[1];
sx q[1];
rz(-1.5473066) q[1];
rz(1.5958105) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(2.3576095) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2326956) q[0];
sx q[0];
rz(-0.43982738) q[0];
sx q[0];
rz(0.16172414) q[0];
rz(-3.0555658) q[2];
sx q[2];
rz(-1.9326397) q[2];
sx q[2];
rz(1.3081683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3938013) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(1.5138813) q[1];
x q[2];
rz(-0.46320398) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(-0.9334329) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(1.5110925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33414267) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(-2.6956218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0837299) q[1];
sx q[1];
rz(-1.5340367) q[1];
sx q[1];
rz(-0.68395331) q[1];
rz(-pi) q[2];
rz(-1.9282354) q[3];
sx q[3];
rz(-2.1921428) q[3];
sx q[3];
rz(0.77222094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0788706) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659347) q[0];
sx q[0];
rz(-1.6871534) q[0];
sx q[0];
rz(-0.12136202) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76639558) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(3.0255084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93563423) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(2.2531855) q[1];
x q[2];
rz(-0.58335431) q[3];
sx q[3];
rz(-1.5837985) q[3];
sx q[3];
rz(1.0540203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.3674659) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-2.6828962) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(1.6470439) q[2];
rz(0.17584569) q[3];
sx q[3];
rz(-2.1538215) q[3];
sx q[3];
rz(0.62721696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];