OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(-2.3119976) q[0];
sx q[0];
rz(-0.15396804) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6217123) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(1.8983311) q[0];
rz(0.36891919) q[2];
sx q[2];
rz(-2.2152165) q[2];
sx q[2];
rz(1.8298139) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4664647) q[1];
sx q[1];
rz(-2.8505241) q[1];
sx q[1];
rz(2.9558099) q[1];
rz(-0.015720856) q[3];
sx q[3];
rz(-2.0831046) q[3];
sx q[3];
rz(0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.1738698) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(-0.57463542) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(-1.8992791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984098) q[0];
sx q[0];
rz(-1.3230723) q[0];
sx q[0];
rz(1.7032911) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.785948) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(-1.7807963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6537135) q[1];
sx q[1];
rz(-0.55760819) q[1];
sx q[1];
rz(-0.30456581) q[1];
rz(-0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-0.58369613) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(2.1247991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4749878) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(-0.41390093) q[0];
rz(-pi) q[1];
rz(-0.74480199) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(1.1643861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6779855) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(2.3198818) q[1];
rz(-pi) q[2];
rz(-1.8524283) q[3];
sx q[3];
rz(-1.6179807) q[3];
sx q[3];
rz(-0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.6695492) q[0];
rz(-0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(2.8947815) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26537672) q[0];
sx q[0];
rz(-0.88337684) q[0];
sx q[0];
rz(2.9486297) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3601801) q[2];
sx q[2];
rz(-1.1014551) q[2];
sx q[2];
rz(-0.59553713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88627316) q[1];
sx q[1];
rz(-1.5181932) q[1];
sx q[1];
rz(-0.37087755) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1770583) q[3];
sx q[3];
rz(-2.4871832) q[3];
sx q[3];
rz(0.27458336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43863338) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(-0.066141733) q[0];
x q[1];
rz(0.18647285) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(-1.965167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0548965) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(0.30026786) q[1];
x q[2];
rz(0.31110839) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(-2.6274519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(0.38875368) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(0.1246917) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9799177) q[0];
sx q[0];
rz(-1.3945762) q[0];
sx q[0];
rz(-1.6582703) q[0];
x q[1];
rz(2.7331946) q[2];
sx q[2];
rz(-2.4917267) q[2];
sx q[2];
rz(1.0330531) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2248968) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(-0.7373666) q[1];
rz(-pi) q[2];
rz(2.9051404) q[3];
sx q[3];
rz(-1.7835622) q[3];
sx q[3];
rz(-0.80602431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(2.0986309) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(1.0707062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8541504) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(2.6575412) q[0];
rz(-2.7878739) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(-3.0596717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6003905) q[1];
sx q[1];
rz(-1.6704847) q[1];
sx q[1];
rz(2.3669835) q[1];
rz(1.5504863) q[3];
sx q[3];
rz(-1.2386525) q[3];
sx q[3];
rz(0.16196812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(-1.2063684) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763801) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(3.0239848) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1768441) q[2];
sx q[2];
rz(-1.7272005) q[2];
sx q[2];
rz(2.1795189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7906588) q[1];
sx q[1];
rz(-0.62218636) q[1];
sx q[1];
rz(1.15637) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3341228) q[3];
sx q[3];
rz(-1.8261357) q[3];
sx q[3];
rz(-2.1047999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751574) q[0];
sx q[0];
rz(-0.046910722) q[0];
sx q[0];
rz(-2.5527918) q[0];
rz(-0.75195306) q[2];
sx q[2];
rz(-0.30138902) q[2];
sx q[2];
rz(1.4434659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(1.4700252) q[1];
rz(-0.30236249) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(0.16929786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(2.8607686) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40598665) q[0];
sx q[0];
rz(-1.9970903) q[0];
sx q[0];
rz(-2.6398525) q[0];
rz(-pi) q[1];
rz(-1.7634723) q[2];
sx q[2];
rz(-1.3609876) q[2];
sx q[2];
rz(-0.10373058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22051375) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(-1.6284579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9032352) q[3];
sx q[3];
rz(-1.7566163) q[3];
sx q[3];
rz(-2.9591054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(-2.279083) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
