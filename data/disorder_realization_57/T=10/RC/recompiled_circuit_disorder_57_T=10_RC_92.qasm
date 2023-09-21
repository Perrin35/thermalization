OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(-0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42703585) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(-1.3546076) q[0];
rz(-2.4389624) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(-1.4544912) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4431475) q[1];
sx q[1];
rz(-1.5138211) q[1];
sx q[1];
rz(-0.2351825) q[1];
rz(-0.29294149) q[3];
sx q[3];
rz(-1.5353234) q[3];
sx q[3];
rz(-0.5637416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.6137326) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(-2.3221827) q[0];
rz(-2.8557414) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768515) q[0];
sx q[0];
rz(-0.74155945) q[0];
sx q[0];
rz(-2.5624647) q[0];
x q[1];
rz(1.9752713) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(1.5420367) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5310865) q[1];
sx q[1];
rz(-2.5889581) q[1];
sx q[1];
rz(-2.3081739) q[1];
rz(0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(0.30360046) q[0];
rz(1.7595694) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(2.4940925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30595368) q[0];
sx q[0];
rz(-1.6157955) q[0];
sx q[0];
rz(0.63458981) q[0];
x q[1];
rz(0.60909231) q[2];
sx q[2];
rz(-2.3272772) q[2];
sx q[2];
rz(2.7012205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69933575) q[1];
sx q[1];
rz(-1.934821) q[1];
sx q[1];
rz(1.0815716) q[1];
rz(0.59433523) q[3];
sx q[3];
rz(-2.1777275) q[3];
sx q[3];
rz(2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(-0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340387) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(0.87189829) q[0];
rz(-0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(-0.29104582) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8431393) q[0];
sx q[0];
rz(-1.342134) q[0];
sx q[0];
rz(2.162621) q[0];
rz(-pi) q[1];
rz(-2.8486738) q[2];
sx q[2];
rz(-1.8790763) q[2];
sx q[2];
rz(-0.76383797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93451553) q[1];
sx q[1];
rz(-2.0743437) q[1];
sx q[1];
rz(2.3758923) q[1];
rz(-pi) q[2];
rz(-0.16104161) q[3];
sx q[3];
rz(-2.1412009) q[3];
sx q[3];
rz(-2.963484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3691833) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13957025) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(0.64506662) q[0];
rz(2.5158665) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-2.5114139) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3483) q[0];
sx q[0];
rz(-0.59069809) q[0];
sx q[0];
rz(-2.9212055) q[0];
rz(-0.39406392) q[2];
sx q[2];
rz(-0.85959083) q[2];
sx q[2];
rz(1.994339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9471776) q[1];
sx q[1];
rz(-1.0395323) q[1];
sx q[1];
rz(0.26015014) q[1];
rz(-2.9003521) q[3];
sx q[3];
rz(-2.5706228) q[3];
sx q[3];
rz(-0.2824479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(2.3278918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57665598) q[0];
sx q[0];
rz(-1.7840958) q[0];
sx q[0];
rz(3.0773786) q[0];
rz(-pi) q[1];
rz(-2.5021624) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(0.30314988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4390347) q[1];
sx q[1];
rz(-0.82779373) q[1];
sx q[1];
rz(-2.031416) q[1];
x q[2];
rz(0.86088647) q[3];
sx q[3];
rz(-2.6885899) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-0.6860835) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(-0.66326052) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903117) q[0];
sx q[0];
rz(-0.65070063) q[0];
sx q[0];
rz(1.531321) q[0];
rz(-2.8510423) q[2];
sx q[2];
rz(-2.4834589) q[2];
sx q[2];
rz(1.5302637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2901193) q[1];
sx q[1];
rz(-1.0901325) q[1];
sx q[1];
rz(2.9139247) q[1];
rz(-0.56179629) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(2.3464349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.474581) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(0.015080301) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(-0.10738871) q[0];
rz(-2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(2.6838141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900507) q[0];
sx q[0];
rz(-2.4180782) q[0];
sx q[0];
rz(-1.0735805) q[0];
rz(-pi) q[1];
rz(-2.1019943) q[2];
sx q[2];
rz(-2.2441412) q[2];
sx q[2];
rz(1.4199867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.01955186) q[1];
sx q[1];
rz(-1.8035839) q[1];
sx q[1];
rz(1.9368534) q[1];
x q[2];
rz(-1.9321953) q[3];
sx q[3];
rz(-2.0933588) q[3];
sx q[3];
rz(-0.88025974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(-1.7163904) q[3];
sx q[3];
rz(-2.5148354) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872221) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(-1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(2.5568533) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7493233) q[0];
sx q[0];
rz(-1.682878) q[0];
sx q[0];
rz(-2.6507676) q[0];
x q[1];
rz(-1.936475) q[2];
sx q[2];
rz(-2.0846016) q[2];
sx q[2];
rz(0.50525451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19496275) q[1];
sx q[1];
rz(-1.1319036) q[1];
sx q[1];
rz(-1.4942601) q[1];
rz(-pi) q[2];
rz(-0.14296602) q[3];
sx q[3];
rz(-1.3060095) q[3];
sx q[3];
rz(-0.20204443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(-2.8653223) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(-0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(2.4861091) q[0];
rz(-1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(-1.5225333) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1725384) q[0];
sx q[0];
rz(-1.8803055) q[0];
sx q[0];
rz(-2.4597416) q[0];
rz(-pi) q[1];
rz(-2.4535975) q[2];
sx q[2];
rz(-1.8383467) q[2];
sx q[2];
rz(-0.5984532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22776991) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(0.63197244) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7195815) q[3];
sx q[3];
rz(-1.2000298) q[3];
sx q[3];
rz(1.5164204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15554252) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(-1.6048253) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(2.4189667) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(1.7799829) q[2];
sx q[2];
rz(-1.5968512) q[2];
sx q[2];
rz(-3.130393) q[2];
rz(1.4747254) q[3];
sx q[3];
rz(-0.69312743) q[3];
sx q[3];
rz(3.0467924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];