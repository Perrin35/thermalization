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
rz(2.6101987) q[0];
rz(-2.9040789) q[1];
sx q[1];
rz(-1.3637435) q[1];
sx q[1];
rz(-1.9385424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2837613) q[0];
sx q[0];
rz(-2.8537822) q[0];
sx q[0];
rz(-0.83588155) q[0];
rz(-pi) q[1];
rz(-2.4389624) q[2];
sx q[2];
rz(-0.32877562) q[2];
sx q[2];
rz(1.4544912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0002999) q[1];
sx q[1];
rz(-1.3360026) q[1];
sx q[1];
rz(-1.512212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5337464) q[3];
sx q[3];
rz(-1.2780446) q[3];
sx q[3];
rz(-1.0177514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-2.9425088) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(-2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99877015) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(-2.8557414) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.8751289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931912) q[0];
sx q[0];
rz(-1.9494434) q[0];
sx q[0];
rz(0.65403954) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1663213) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(1.599556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30217583) q[1];
sx q[1];
rz(-1.9315047) q[1];
sx q[1];
rz(1.9990666) q[1];
rz(-pi) q[2];
rz(-0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(-1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(-0.087163838) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.7595694) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(2.4940925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.203813) q[0];
sx q[0];
rz(-2.505629) q[0];
sx q[0];
rz(3.0657835) q[0];
rz(-pi) q[1];
rz(-2.5325003) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(-2.7012205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4574617) q[1];
sx q[1];
rz(-2.0254454) q[1];
sx q[1];
rz(-2.734114) q[1];
rz(-pi) q[2];
rz(-0.59433523) q[3];
sx q[3];
rz(-0.96386516) q[3];
sx q[3];
rz(2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5807318) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(-0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(0.87189829) q[0];
rz(2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(-2.8505468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211348) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(0.27340425) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8918858) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(2.2433787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17109045) q[1];
sx q[1];
rz(-2.2541752) q[1];
sx q[1];
rz(0.6716397) q[1];
rz(1.8157186) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(0.11412379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(2.496526) q[0];
rz(-2.5158665) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-0.63017875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1030578) q[0];
sx q[0];
rz(-1.6928506) q[0];
sx q[0];
rz(-0.5794258) q[0];
x q[1];
rz(-0.39406392) q[2];
sx q[2];
rz(-2.2820018) q[2];
sx q[2];
rz(1.1472536) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2423357) q[1];
sx q[1];
rz(-1.3471654) q[1];
sx q[1];
rz(-1.0244589) q[1];
x q[2];
rz(1.7230677) q[3];
sx q[3];
rz(-2.1232743) q[3];
sx q[3];
rz(-0.0020364062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(-0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(0.71682799) q[0];
rz(-0.43771276) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(-2.3278918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5649367) q[0];
sx q[0];
rz(-1.7840958) q[0];
sx q[0];
rz(-0.06421406) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4196163) q[2];
sx q[2];
rz(-2.0275896) q[2];
sx q[2];
rz(1.4095838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5972283) q[1];
sx q[1];
rz(-1.2372984) q[1];
sx q[1];
rz(-2.3436105) q[1];
rz(-pi) q[2];
rz(0.86088647) q[3];
sx q[3];
rz(-2.6885899) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.1876594) q[2];
rz(2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903117) q[0];
sx q[0];
rz(-0.65070063) q[0];
sx q[0];
rz(-1.6102717) q[0];
rz(-1.7887605) q[2];
sx q[2];
rz(-0.94467615) q[2];
sx q[2];
rz(1.2499714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85147334) q[1];
sx q[1];
rz(-1.0901325) q[1];
sx q[1];
rz(0.22766797) q[1];
rz(-2.0459941) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(-1.0218395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(-0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(3.0342039) q[0];
rz(0.30474162) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-0.4577786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37516884) q[0];
sx q[0];
rz(-0.94978118) q[0];
sx q[0];
rz(0.39874886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7465676) q[2];
sx q[2];
rz(-1.1636359) q[2];
sx q[2];
rz(0.50200576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0929872) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(2.155817) q[1];
rz(2.5896766) q[3];
sx q[3];
rz(-1.2593927) q[3];
sx q[3];
rz(0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.8553998) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872221) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(-1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-0.58473933) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9033602) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(1.4438629) q[0];
x q[1];
rz(-1.2051177) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(0.50525451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4084088) q[1];
sx q[1];
rz(-1.5015263) q[1];
sx q[1];
rz(2.7015711) q[1];
rz(-pi) q[2];
rz(-0.14296602) q[3];
sx q[3];
rz(-1.8355832) q[3];
sx q[3];
rz(0.20204443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(-0.27627036) q[2];
rz(-2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(-2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722512) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(-1.6190593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(1.9615016) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40752878) q[2];
sx q[2];
rz(-2.4113843) q[2];
sx q[2];
rz(-0.66116316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9138227) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(-0.63197244) q[1];
rz(-pi) q[2];
rz(-1.9736245) q[3];
sx q[3];
rz(-1.1790923) q[3];
sx q[3];
rz(-0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15554252) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(2.299451) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(0.026528927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-0.72262598) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(-1.6956383) q[2];
sx q[2];
rz(-0.21077934) q[2];
sx q[2];
rz(-1.6817033) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];