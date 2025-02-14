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
rz(1.0951618) q[0];
sx q[0];
rz(0.14552966) q[0];
sx q[0];
rz(9.9280823) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(-1.6961179) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2615525) q[0];
sx q[0];
rz(-1.6140811) q[0];
sx q[0];
rz(-2.7425984) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93930556) q[2];
sx q[2];
rz(-0.92512265) q[2];
sx q[2];
rz(1.948057) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9933982) q[1];
sx q[1];
rz(-0.96635093) q[1];
sx q[1];
rz(1.1682214) q[1];
x q[2];
rz(1.5863933) q[3];
sx q[3];
rz(-1.6199779) q[3];
sx q[3];
rz(1.5893857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4989) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(-1.6373681) q[2];
rz(2.4046992) q[3];
sx q[3];
rz(-1.5259909) q[3];
sx q[3];
rz(-0.98114291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7351643) q[0];
sx q[0];
rz(-2.6728215) q[0];
sx q[0];
rz(-2.6306187) q[0];
rz(-1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-2.8151292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0827507) q[0];
sx q[0];
rz(-1.8371757) q[0];
sx q[0];
rz(-2.8681115) q[0];
rz(-pi) q[1];
rz(2.9316177) q[2];
sx q[2];
rz(-0.89824235) q[2];
sx q[2];
rz(2.5598516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6447811) q[1];
sx q[1];
rz(-2.4678951) q[1];
sx q[1];
rz(1.0717057) q[1];
rz(2.191698) q[3];
sx q[3];
rz(-2.6767459) q[3];
sx q[3];
rz(-0.45528938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9713126) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(1.2962606) q[2];
rz(-0.41804677) q[3];
sx q[3];
rz(-2.1635735) q[3];
sx q[3];
rz(0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0728077) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(0.54247722) q[0];
rz(1.7659278) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-0.03104041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.870793) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(0.55678456) q[0];
rz(-2.9152053) q[2];
sx q[2];
rz(-2.2157266) q[2];
sx q[2];
rz(-1.8075391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2636288) q[1];
sx q[1];
rz(-2.1641556) q[1];
sx q[1];
rz(-2.8441162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5543082) q[3];
sx q[3];
rz(-1.015839) q[3];
sx q[3];
rz(-2.6944427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5151908) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(2.314563) q[2];
rz(2.6229897) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(-1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14722918) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.2116785) q[0];
rz(-1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.3406219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3723243) q[0];
sx q[0];
rz(-0.71592531) q[0];
sx q[0];
rz(0.47874079) q[0];
rz(-pi) q[1];
rz(-1.421809) q[2];
sx q[2];
rz(-1.8049311) q[2];
sx q[2];
rz(0.17727597) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27599469) q[1];
sx q[1];
rz(-1.0038687) q[1];
sx q[1];
rz(-1.2204942) q[1];
rz(-2.2068437) q[3];
sx q[3];
rz(-0.49643653) q[3];
sx q[3];
rz(1.2477341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2010605) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-0.96735668) q[3];
sx q[3];
rz(0.78260666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56213266) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(2.539047) q[0];
rz(-2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(-0.7811195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21602042) q[0];
sx q[0];
rz(-1.7402152) q[0];
sx q[0];
rz(-1.8623307) q[0];
rz(-pi) q[1];
rz(0.76483043) q[2];
sx q[2];
rz(-2.0553608) q[2];
sx q[2];
rz(0.35540211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.994532) q[1];
sx q[1];
rz(-2.414738) q[1];
sx q[1];
rz(1.4406073) q[1];
rz(-0.1991462) q[3];
sx q[3];
rz(-2.4059787) q[3];
sx q[3];
rz(2.3922331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96131229) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(-0.22892924) q[2];
rz(-1.3198352) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(0.84407097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9127386) q[0];
sx q[0];
rz(-1.5541394) q[0];
sx q[0];
rz(1.9629021) q[0];
rz(-1.9121869) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(-1.8454525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73979171) q[0];
sx q[0];
rz(-1.8515046) q[0];
sx q[0];
rz(-0.14174353) q[0];
rz(-pi) q[1];
rz(0.8528233) q[2];
sx q[2];
rz(-2.1666424) q[2];
sx q[2];
rz(-0.15534523) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78627045) q[1];
sx q[1];
rz(-0.67364022) q[1];
sx q[1];
rz(-2.4060529) q[1];
x q[2];
rz(0.48719897) q[3];
sx q[3];
rz(-1.6604843) q[3];
sx q[3];
rz(-1.439217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65427762) q[2];
sx q[2];
rz(-2.2389905) q[2];
sx q[2];
rz(0.26228341) q[2];
rz(-0.30019635) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(-2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4406141) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(1.3025008) q[0];
rz(-2.473096) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(-2.1065333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651606) q[0];
sx q[0];
rz(-0.4598099) q[0];
sx q[0];
rz(2.0681429) q[0];
rz(-1.9403753) q[2];
sx q[2];
rz(-1.3962922) q[2];
sx q[2];
rz(0.39350393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25495249) q[1];
sx q[1];
rz(-1.4952085) q[1];
sx q[1];
rz(-1.9151389) q[1];
x q[2];
rz(-2.1627102) q[3];
sx q[3];
rz(-0.91875263) q[3];
sx q[3];
rz(2.121821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11789007) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(1.6560076) q[2];
rz(-0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(2.7070467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(-0.27454141) q[0];
rz(1.2046332) q[1];
sx q[1];
rz(-2.5614673) q[1];
sx q[1];
rz(-2.5801632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12073853) q[0];
sx q[0];
rz(-2.0971813) q[0];
sx q[0];
rz(0.83197439) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8738633) q[2];
sx q[2];
rz(-2.2243739) q[2];
sx q[2];
rz(0.67959626) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2711948) q[1];
sx q[1];
rz(-1.8742838) q[1];
sx q[1];
rz(1.1417901) q[1];
rz(-pi) q[2];
rz(-1.2227603) q[3];
sx q[3];
rz(-1.3216747) q[3];
sx q[3];
rz(0.083569738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56208912) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(1.7378463) q[2];
rz(-1.7956519) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(0.48817202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.76029921) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(-2.2316933) q[0];
rz(-1.9445885) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(-0.88814703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67535454) q[0];
sx q[0];
rz(-0.14396891) q[0];
sx q[0];
rz(0.86021249) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0641563) q[2];
sx q[2];
rz(-0.66129035) q[2];
sx q[2];
rz(-1.3232302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5822118) q[1];
sx q[1];
rz(-2.191698) q[1];
sx q[1];
rz(-1.0889174) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4481285) q[3];
sx q[3];
rz(-2.8438206) q[3];
sx q[3];
rz(0.55469251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1278648) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(1.0825276) q[2];
rz(2.9108289) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.3103264) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(0.58151522) q[0];
rz(0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(1.2967348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0688419) q[0];
sx q[0];
rz(-1.5889935) q[0];
sx q[0];
rz(-0.010404603) q[0];
rz(-0.34411547) q[2];
sx q[2];
rz(-1.5617579) q[2];
sx q[2];
rz(-1.2492385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77380005) q[1];
sx q[1];
rz(-1.3145216) q[1];
sx q[1];
rz(-1.132375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6526645) q[3];
sx q[3];
rz(-0.73344123) q[3];
sx q[3];
rz(-1.5239925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41946188) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(-2.6922743) q[2];
rz(-2.8214473) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768059) q[0];
sx q[0];
rz(-2.4644485) q[0];
sx q[0];
rz(1.2039536) q[0];
rz(1.3846579) q[1];
sx q[1];
rz(-0.23778267) q[1];
sx q[1];
rz(2.3429088) q[1];
rz(0.25945406) q[2];
sx q[2];
rz(-2.1726407) q[2];
sx q[2];
rz(3.0214027) q[2];
rz(2.7560227) q[3];
sx q[3];
rz(-1.1299993) q[3];
sx q[3];
rz(-0.53408505) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
