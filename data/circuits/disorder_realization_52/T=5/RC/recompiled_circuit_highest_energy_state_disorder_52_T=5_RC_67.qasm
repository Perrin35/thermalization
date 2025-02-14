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
rz(-1.4327383) q[0];
sx q[0];
rz(-0.32810768) q[0];
sx q[0];
rz(-2.3018667) q[0];
rz(-1.7307164) q[1];
sx q[1];
rz(-1.6003992) q[1];
sx q[1];
rz(2.3720001) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107723) q[0];
sx q[0];
rz(-1.4913173) q[0];
sx q[0];
rz(-3.1067418) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5430449) q[2];
sx q[2];
rz(-0.74615462) q[2];
sx q[2];
rz(-1.0026102) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.283764) q[1];
sx q[1];
rz(-0.35175465) q[1];
sx q[1];
rz(1.5689538) q[1];
rz(2.493164) q[3];
sx q[3];
rz(-2.3852013) q[3];
sx q[3];
rz(-1.7738916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8771693) q[2];
sx q[2];
rz(-1.8921655) q[2];
sx q[2];
rz(-1.9682311) q[2];
rz(0.42826432) q[3];
sx q[3];
rz(-0.56848017) q[3];
sx q[3];
rz(-0.65304023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1332755) q[0];
sx q[0];
rz(-0.44047099) q[0];
sx q[0];
rz(2.0793197) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(1.0947469) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0355201) q[0];
sx q[0];
rz(-0.2600593) q[0];
sx q[0];
rz(-1.3706657) q[0];
x q[1];
rz(0.89265426) q[2];
sx q[2];
rz(-1.2142748) q[2];
sx q[2];
rz(-1.4068436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96910406) q[1];
sx q[1];
rz(-1.8575694) q[1];
sx q[1];
rz(0.4003093) q[1];
rz(-1.7844438) q[3];
sx q[3];
rz(-2.3348138) q[3];
sx q[3];
rz(2.2329604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67794472) q[2];
sx q[2];
rz(-0.28767583) q[2];
sx q[2];
rz(-2.2410683) q[2];
rz(2.1053704) q[3];
sx q[3];
rz(-0.13768727) q[3];
sx q[3];
rz(-3.1305967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66505945) q[0];
sx q[0];
rz(-1.1256951) q[0];
sx q[0];
rz(-1.6735459) q[0];
rz(0.92848575) q[1];
sx q[1];
rz(-1.0037582) q[1];
sx q[1];
rz(-1.4220062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8385561) q[0];
sx q[0];
rz(-1.3232797) q[0];
sx q[0];
rz(-1.4776138) q[0];
x q[1];
rz(2.3238238) q[2];
sx q[2];
rz(-1.0330794) q[2];
sx q[2];
rz(1.0420711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8633008) q[1];
sx q[1];
rz(-1.2108786) q[1];
sx q[1];
rz(0.57508303) q[1];
rz(-pi) q[2];
rz(-0.57392759) q[3];
sx q[3];
rz(-2.6477154) q[3];
sx q[3];
rz(0.18791325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15098393) q[2];
sx q[2];
rz(-1.7093806) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(-3.0153583) q[3];
sx q[3];
rz(-1.9641967) q[3];
sx q[3];
rz(1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6821297) q[0];
sx q[0];
rz(-1.4012902) q[0];
sx q[0];
rz(1.6023741) q[0];
rz(1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.1702671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4919969) q[0];
sx q[0];
rz(-1.2536712) q[0];
sx q[0];
rz(-1.9673504) q[0];
x q[1];
rz(-0.2601461) q[2];
sx q[2];
rz(-2.2739774) q[2];
sx q[2];
rz(1.2728572) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6017157) q[1];
sx q[1];
rz(-1.6605494) q[1];
sx q[1];
rz(1.6513157) q[1];
rz(-pi) q[2];
rz(-0.37335124) q[3];
sx q[3];
rz(-1.7718833) q[3];
sx q[3];
rz(1.1652077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22152659) q[2];
sx q[2];
rz(-2.2130794) q[2];
sx q[2];
rz(-3.0583734) q[2];
rz(-2.4890238) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(2.2799802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457152) q[0];
sx q[0];
rz(-0.28384122) q[0];
sx q[0];
rz(0.19677095) q[0];
rz(-1.1605877) q[1];
sx q[1];
rz(-0.44191688) q[1];
sx q[1];
rz(1.7534076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3144238) q[0];
sx q[0];
rz(-1.4008796) q[0];
sx q[0];
rz(-1.569237) q[0];
rz(-pi) q[1];
rz(0.20241995) q[2];
sx q[2];
rz(-1.4223863) q[2];
sx q[2];
rz(-2.9925516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3922351) q[1];
sx q[1];
rz(-0.92441407) q[1];
sx q[1];
rz(-2.9254167) q[1];
x q[2];
rz(-0.36404534) q[3];
sx q[3];
rz(-0.87556404) q[3];
sx q[3];
rz(-2.1882868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5768726) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(-1.2811071) q[2];
rz(2.7992904) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(0.92882338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-1.0973278) q[0];
sx q[0];
rz(2.1088364) q[0];
rz(-2.4646941) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(-2.3597609) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9932258) q[0];
sx q[0];
rz(-1.2162104) q[0];
sx q[0];
rz(-0.605159) q[0];
rz(-pi) q[1];
x q[1];
rz(1.621945) q[2];
sx q[2];
rz(-1.5384073) q[2];
sx q[2];
rz(0.50078228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1429174) q[1];
sx q[1];
rz(-1.8774722) q[1];
sx q[1];
rz(-2.5481413) q[1];
x q[2];
rz(0.12198066) q[3];
sx q[3];
rz(-2.8013419) q[3];
sx q[3];
rz(-1.550479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8402164) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(2.2678579) q[2];
rz(0.78911632) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9588722) q[0];
sx q[0];
rz(-3.0954269) q[0];
sx q[0];
rz(0.090855457) q[0];
rz(2.5317522) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-2.5441433) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0598568) q[0];
sx q[0];
rz(-1.6993521) q[0];
sx q[0];
rz(0.14692326) q[0];
rz(-1.1979073) q[2];
sx q[2];
rz(-1.2593566) q[2];
sx q[2];
rz(2.5562037) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8272463) q[1];
sx q[1];
rz(-1.1097849) q[1];
sx q[1];
rz(-0.48996144) q[1];
rz(3.1037056) q[3];
sx q[3];
rz(-1.1514613) q[3];
sx q[3];
rz(-1.9333206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7647543) q[2];
sx q[2];
rz(-2.59616) q[2];
sx q[2];
rz(0.1405912) q[2];
rz(1.2815255) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-0.52711058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22290467) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(1.5284398) q[0];
rz(-2.3279922) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(-2.334107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59283108) q[0];
sx q[0];
rz(-0.27540311) q[0];
sx q[0];
rz(-2.3345678) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3027655) q[2];
sx q[2];
rz(-0.45739528) q[2];
sx q[2];
rz(-3.0283324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1450913) q[1];
sx q[1];
rz(-1.1942625) q[1];
sx q[1];
rz(2.242617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0309615) q[3];
sx q[3];
rz(-1.9412517) q[3];
sx q[3];
rz(-2.3453804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7878824) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(0.14459571) q[2];
rz(-1.1041798) q[3];
sx q[3];
rz(-0.2140597) q[3];
sx q[3];
rz(-2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5398194) q[0];
sx q[0];
rz(-1.1310534) q[0];
sx q[0];
rz(2.5850776) q[0];
rz(-0.76983184) q[1];
sx q[1];
rz(-3.0172805) q[1];
sx q[1];
rz(2.6803023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3614013) q[0];
sx q[0];
rz(-2.5690329) q[0];
sx q[0];
rz(1.6587692) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3084505) q[2];
sx q[2];
rz(-2.1066446) q[2];
sx q[2];
rz(1.3539202) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76595054) q[1];
sx q[1];
rz(-1.2772755) q[1];
sx q[1];
rz(2.3883345) q[1];
rz(-2.5784745) q[3];
sx q[3];
rz(-0.65059911) q[3];
sx q[3];
rz(1.4103149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43872675) q[2];
sx q[2];
rz(-0.30539572) q[2];
sx q[2];
rz(-0.76496441) q[2];
rz(1.2139828) q[3];
sx q[3];
rz(-2.5524804) q[3];
sx q[3];
rz(0.37863076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42534378) q[0];
sx q[0];
rz(-3.0916164) q[0];
sx q[0];
rz(2.7783527) q[0];
rz(-1.7323469) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(-0.22663103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55174819) q[0];
sx q[0];
rz(-1.4355816) q[0];
sx q[0];
rz(1.5677794) q[0];
rz(-0.93904597) q[2];
sx q[2];
rz(-2.2151196) q[2];
sx q[2];
rz(-2.7959888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5780826) q[1];
sx q[1];
rz(-0.77111608) q[1];
sx q[1];
rz(-2.0457436) q[1];
rz(3.111242) q[3];
sx q[3];
rz(-0.33799809) q[3];
sx q[3];
rz(-2.3299467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54214415) q[2];
sx q[2];
rz(-2.6207974) q[2];
sx q[2];
rz(-2.4929872) q[2];
rz(-3.0829698) q[3];
sx q[3];
rz(-2.5128745) q[3];
sx q[3];
rz(0.64529836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36515737) q[0];
sx q[0];
rz(-1.9226274) q[0];
sx q[0];
rz(2.6489039) q[0];
rz(1.0538712) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(-1.8070167) q[2];
sx q[2];
rz(-1.5945572) q[2];
sx q[2];
rz(-0.40785892) q[2];
rz(1.1751997) q[3];
sx q[3];
rz(-2.4064872) q[3];
sx q[3];
rz(-0.088633782) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
