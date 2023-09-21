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
rz(0.33831236) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96903893) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(2.88455) q[0];
rz(2.7726735) q[2];
sx q[2];
rz(-0.92637617) q[2];
sx q[2];
rz(1.8298139) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6602064) q[1];
sx q[1];
rz(-1.8567137) q[1];
sx q[1];
rz(-1.5155161) q[1];
rz(-pi) q[2];
rz(0.015720856) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-3.0766292) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.2423135) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500538) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(0.48135249) q[0];
rz(-pi) q[1];
rz(1.3556446) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(-1.7807963) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9638483) q[1];
sx q[1];
rz(-1.7301534) q[1];
sx q[1];
rz(2.6049155) q[1];
rz(-0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(0.58369613) q[2];
rz(-0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-2.4131391) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-2.1247991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66660488) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(2.7276917) q[0];
rz(-1.0810034) q[2];
sx q[2];
rz(-2.042633) q[2];
sx q[2];
rz(-2.8420198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9079202) q[1];
sx q[1];
rz(-2.3111812) q[1];
sx q[1];
rz(2.9552712) q[1];
rz(-pi) q[2];
rz(0.049116491) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(-2.1360872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8762159) q[0];
sx q[0];
rz(-0.88337684) q[0];
sx q[0];
rz(2.9486297) q[0];
rz(-2.6631782) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(-1.071655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66407953) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(-1.627229) q[1];
rz(-pi) q[2];
rz(2.187192) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43863338) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(3.0754509) q[0];
rz(0.18647285) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(1.1764256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3779113) q[1];
sx q[1];
rz(-1.3994201) q[1];
sx q[1];
rz(2.548449) q[1];
rz(-pi) q[2];
rz(-1.3316657) q[3];
sx q[3];
rz(-0.93591792) q[3];
sx q[3];
rz(0.90415108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4437618) q[0];
sx q[0];
rz(-0.19653453) q[0];
sx q[0];
rz(-2.6854808) q[0];
rz(2.5325534) q[2];
sx q[2];
rz(-1.3281203) q[2];
sx q[2];
rz(-2.2720624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71619892) q[1];
sx q[1];
rz(-2.4024706) q[1];
sx q[1];
rz(3.0576586) q[1];
rz(-2.3966339) q[3];
sx q[3];
rz(-2.8248441) q[3];
sx q[3];
rz(-1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-2.0708864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2830414) q[0];
sx q[0];
rz(-2.0548477) q[0];
sx q[0];
rz(-1.5714684) q[0];
x q[1];
rz(-2.7942065) q[2];
sx q[2];
rz(-1.6958478) q[2];
sx q[2];
rz(1.8206247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.067498265) q[1];
sx q[1];
rz(-0.80103445) q[1];
sx q[1];
rz(-1.7098411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5911063) q[3];
sx q[3];
rz(-1.9029402) q[3];
sx q[3];
rz(0.16196812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(2.1941197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8038717) q[0];
sx q[0];
rz(-0.54715711) q[0];
sx q[0];
rz(-1.3756115) q[0];
rz(1.7296373) q[2];
sx q[2];
rz(-1.3961332) q[2];
sx q[2];
rz(-0.58089248) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87654963) q[1];
sx q[1];
rz(-1.8076841) q[1];
sx q[1];
rz(0.9898647) q[1];
x q[2];
rz(-0.8074699) q[3];
sx q[3];
rz(-1.8261357) q[3];
sx q[3];
rz(-2.1047999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(0.20283094) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3644667) q[0];
sx q[0];
rz(-1.5317894) q[0];
sx q[0];
rz(1.5447306) q[0];
rz(-pi) q[1];
rz(-2.9183396) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(-0.60165652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(-1.6715675) q[1];
x q[2];
rz(-0.30236249) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.330864) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(-0.75673639) q[0];
rz(-pi) q[1];
rz(0.73239399) q[2];
sx q[2];
rz(-0.28389441) q[2];
sx q[2];
rz(-2.4925799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.909543) q[1];
sx q[1];
rz(-0.068040158) q[1];
sx q[1];
rz(-1.0104936) q[1];
rz(-1.0481846) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(-2.2446333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(0.87396809) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(0.59861029) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(-0.35691805) q[3];
sx q[3];
rz(-1.1834984) q[3];
sx q[3];
rz(-0.055565861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];