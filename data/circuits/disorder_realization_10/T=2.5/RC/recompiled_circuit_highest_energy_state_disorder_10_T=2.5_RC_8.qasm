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
rz(2.4438357) q[0];
sx q[0];
rz(-1.9480167) q[0];
sx q[0];
rz(0.20456631) q[0];
rz(2.3808631) q[1];
sx q[1];
rz(-1.7525571) q[1];
sx q[1];
rz(-1.3995481) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4142864) q[0];
sx q[0];
rz(-0.23972962) q[0];
sx q[0];
rz(-1.1778465) q[0];
rz(-pi) q[1];
rz(2.7367107) q[2];
sx q[2];
rz(-1.153115) q[2];
sx q[2];
rz(1.6987125) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5198361) q[1];
sx q[1];
rz(-0.43070983) q[1];
sx q[1];
rz(-2.9309896) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70943009) q[3];
sx q[3];
rz(-1.3660079) q[3];
sx q[3];
rz(3.0258609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(0.48933634) q[2];
rz(-2.1183744) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(-2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045687549) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(3.0838222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9607769) q[0];
sx q[0];
rz(-1.0313927) q[0];
sx q[0];
rz(-2.8404854) q[0];
rz(-pi) q[1];
rz(0.34681706) q[2];
sx q[2];
rz(-2.8379734) q[2];
sx q[2];
rz(-0.70250073) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7747468) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(-2.0319035) q[1];
x q[2];
rz(-2.9936419) q[3];
sx q[3];
rz(-2.2770364) q[3];
sx q[3];
rz(2.1220292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94566655) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(1.2954953) q[2];
rz(-0.96674353) q[3];
sx q[3];
rz(-2.3714378) q[3];
sx q[3];
rz(-2.3841592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029723786) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(-1.7080074) q[0];
rz(-3.0729821) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(-2.5624018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5369956) q[0];
sx q[0];
rz(-1.4791577) q[0];
sx q[0];
rz(3.0834404) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0036598031) q[2];
sx q[2];
rz(-2.6658667) q[2];
sx q[2];
rz(1.9329247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3855558) q[1];
sx q[1];
rz(-1.3446523) q[1];
sx q[1];
rz(-1.534034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0924358) q[3];
sx q[3];
rz(-1.1598827) q[3];
sx q[3];
rz(2.876407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8706943) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(2.9581621) q[2];
rz(-2.3373248) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(2.542069) q[0];
rz(0.41246688) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(1.0106769) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090792716) q[0];
sx q[0];
rz(-2.7677892) q[0];
sx q[0];
rz(2.7888377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.050881906) q[2];
sx q[2];
rz(-0.93928277) q[2];
sx q[2];
rz(-0.33755195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4356723) q[1];
sx q[1];
rz(-1.8413891) q[1];
sx q[1];
rz(2.1822004) q[1];
rz(-2.385731) q[3];
sx q[3];
rz(-1.6718277) q[3];
sx q[3];
rz(-1.3974691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7708873) q[2];
sx q[2];
rz(-2.7949896) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(0.62234771) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(-2.5573825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36251003) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(-2.1146178) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(-1.2170353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.352223) q[0];
sx q[0];
rz(-0.93607157) q[0];
sx q[0];
rz(-0.85294368) q[0];
rz(1.4929068) q[2];
sx q[2];
rz(-2.1078034) q[2];
sx q[2];
rz(-0.80918771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0266307) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(-0.86905713) q[1];
rz(1.1547791) q[3];
sx q[3];
rz(-2.3868594) q[3];
sx q[3];
rz(2.5861135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7050742) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-0.77318937) q[2];
rz(3.0298722) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(-2.2022061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(-2.7146085) q[0];
rz(-2.2110979) q[1];
sx q[1];
rz(-3.1246779) q[1];
sx q[1];
rz(-2.6771136) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8278566) q[0];
sx q[0];
rz(-1.858874) q[0];
sx q[0];
rz(-1.7749857) q[0];
x q[1];
rz(2.802425) q[2];
sx q[2];
rz(-0.99771032) q[2];
sx q[2];
rz(3.1011875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.528094) q[1];
sx q[1];
rz(-1.6069176) q[1];
sx q[1];
rz(-0.36904676) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1464981) q[3];
sx q[3];
rz(-1.5992377) q[3];
sx q[3];
rz(2.6928296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6102607) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(0.28826928) q[2];
rz(1.0432976) q[3];
sx q[3];
rz(-2.5259924) q[3];
sx q[3];
rz(-2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6601324) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(0.75755358) q[0];
rz(-3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(0.048197897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.165503) q[0];
sx q[0];
rz(-2.5565326) q[0];
sx q[0];
rz(1.5173605) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2343302) q[2];
sx q[2];
rz(-2.5162553) q[2];
sx q[2];
rz(-0.71695527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0695677) q[1];
sx q[1];
rz(-1.9981019) q[1];
sx q[1];
rz(-0.91520379) q[1];
x q[2];
rz(-0.052560135) q[3];
sx q[3];
rz(-1.1014928) q[3];
sx q[3];
rz(1.8183501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66182071) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(3.0721967) q[2];
rz(-1.5540468) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(-0.21849304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7540392) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(1.4132502) q[0];
rz(0.82855254) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(2.5989596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4601645) q[0];
sx q[0];
rz(-1.6896473) q[0];
sx q[0];
rz(-0.14560856) q[0];
rz(-pi) q[1];
rz(1.2130402) q[2];
sx q[2];
rz(-1.2160436) q[2];
sx q[2];
rz(0.32217978) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5396351) q[1];
sx q[1];
rz(-1.4550721) q[1];
sx q[1];
rz(1.6249815) q[1];
x q[2];
rz(-0.43482921) q[3];
sx q[3];
rz(-1.2246338) q[3];
sx q[3];
rz(2.1157672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98216206) q[2];
sx q[2];
rz(-0.40951481) q[2];
sx q[2];
rz(0.25992599) q[2];
rz(1.0204756) q[3];
sx q[3];
rz(-0.25769886) q[3];
sx q[3];
rz(-0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(1.4790685) q[0];
rz(-1.5326477) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(-2.398568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8296513) q[0];
sx q[0];
rz(-0.52642979) q[0];
sx q[0];
rz(1.3408324) q[0];
rz(0.9833588) q[2];
sx q[2];
rz(-2.5134216) q[2];
sx q[2];
rz(-0.52631015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8240717) q[1];
sx q[1];
rz(-1.5050384) q[1];
sx q[1];
rz(-0.055067896) q[1];
rz(-1.6017655) q[3];
sx q[3];
rz(-1.4659229) q[3];
sx q[3];
rz(-2.7978143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4674025) q[2];
sx q[2];
rz(-2.3154066) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.2286681) q[3];
sx q[3];
rz(-2.3384371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(0.75206494) q[0];
rz(1.1532785) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(-0.28336743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8077791) q[0];
sx q[0];
rz(-0.54610836) q[0];
sx q[0];
rz(2.9696652) q[0];
x q[1];
rz(0.81895639) q[2];
sx q[2];
rz(-1.6791846) q[2];
sx q[2];
rz(-0.10486952) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6436037) q[1];
sx q[1];
rz(-0.97103968) q[1];
sx q[1];
rz(1.1685755) q[1];
rz(-pi) q[2];
rz(1.9612736) q[3];
sx q[3];
rz(-0.9538528) q[3];
sx q[3];
rz(-2.4882567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2470384) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(1.438633) q[2];
rz(-0.29397193) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033584874) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(-1.2106566) q[2];
sx q[2];
rz(-1.4317206) q[2];
sx q[2];
rz(1.7862873) q[2];
rz(2.7593437) q[3];
sx q[3];
rz(-1.8235689) q[3];
sx q[3];
rz(0.27128661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
