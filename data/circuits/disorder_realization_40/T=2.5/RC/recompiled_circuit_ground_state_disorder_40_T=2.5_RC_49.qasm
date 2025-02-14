OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7933529) q[0];
sx q[0];
rz(-1.5433595) q[0];
sx q[0];
rz(-1.6399075) q[0];
rz(-1.3950672) q[1];
sx q[1];
rz(2.7434064) q[1];
sx q[1];
rz(12.413496) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81992309) q[0];
sx q[0];
rz(-1.7060486) q[0];
sx q[0];
rz(3.0514984) q[0];
rz(-pi) q[1];
rz(-0.35043535) q[2];
sx q[2];
rz(-1.7261243) q[2];
sx q[2];
rz(2.2969579) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4248878) q[1];
sx q[1];
rz(-0.198303) q[1];
sx q[1];
rz(1.2601869) q[1];
x q[2];
rz(-2.2431668) q[3];
sx q[3];
rz(-1.200496) q[3];
sx q[3];
rz(2.680925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.07831002) q[2];
sx q[2];
rz(-1.4274884) q[2];
sx q[2];
rz(-0.59986344) q[2];
rz(2.6317224) q[3];
sx q[3];
rz(-0.32250914) q[3];
sx q[3];
rz(2.9738284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55709368) q[0];
sx q[0];
rz(-2.2826513) q[0];
sx q[0];
rz(0.14228819) q[0];
rz(1.2917057) q[1];
sx q[1];
rz(-0.53304356) q[1];
sx q[1];
rz(-0.60089111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639524) q[0];
sx q[0];
rz(-0.48177347) q[0];
sx q[0];
rz(-1.0642306) q[0];
x q[1];
rz(0.49111207) q[2];
sx q[2];
rz(-0.69780707) q[2];
sx q[2];
rz(1.5581824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3004724) q[1];
sx q[1];
rz(-2.4266647) q[1];
sx q[1];
rz(0.27473533) q[1];
x q[2];
rz(1.7606335) q[3];
sx q[3];
rz(-2.077436) q[3];
sx q[3];
rz(1.1200126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94704023) q[2];
sx q[2];
rz(-1.9779454) q[2];
sx q[2];
rz(-2.2770503) q[2];
rz(-2.7583127) q[3];
sx q[3];
rz(-2.7195103) q[3];
sx q[3];
rz(-2.2569136) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8466723) q[0];
sx q[0];
rz(-2.8546794) q[0];
sx q[0];
rz(-0.82557803) q[0];
rz(-2.3041252) q[1];
sx q[1];
rz(-1.0191963) q[1];
sx q[1];
rz(2.3501863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3986555) q[0];
sx q[0];
rz(-0.81327945) q[0];
sx q[0];
rz(-0.11095993) q[0];
rz(-2.2242658) q[2];
sx q[2];
rz(-2.978108) q[2];
sx q[2];
rz(1.2965073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36230642) q[1];
sx q[1];
rz(-2.0842881) q[1];
sx q[1];
rz(-1.1720285) q[1];
rz(-1.0867001) q[3];
sx q[3];
rz(-1.181385) q[3];
sx q[3];
rz(2.2661346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6174263) q[2];
sx q[2];
rz(-2.7238621) q[2];
sx q[2];
rz(1.8903271) q[2];
rz(1.9620126) q[3];
sx q[3];
rz(-1.1359295) q[3];
sx q[3];
rz(2.484926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11611045) q[0];
sx q[0];
rz(-1.9132834) q[0];
sx q[0];
rz(2.3226698) q[0];
rz(-3.0602835) q[1];
sx q[1];
rz(-2.5550877) q[1];
sx q[1];
rz(-0.78048817) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1063447) q[0];
sx q[0];
rz(-0.37899932) q[0];
sx q[0];
rz(2.9826791) q[0];
x q[1];
rz(-2.221406) q[2];
sx q[2];
rz(-1.4942385) q[2];
sx q[2];
rz(-0.45659846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8185212) q[1];
sx q[1];
rz(-1.8667606) q[1];
sx q[1];
rz(0.37533111) q[1];
x q[2];
rz(-1.9742161) q[3];
sx q[3];
rz(-1.2781004) q[3];
sx q[3];
rz(3.0967874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2180194) q[2];
sx q[2];
rz(-1.6417445) q[2];
sx q[2];
rz(1.8531331) q[2];
rz(-1.4675379) q[3];
sx q[3];
rz(-0.54687423) q[3];
sx q[3];
rz(0.43681496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9913637) q[0];
sx q[0];
rz(-0.65058351) q[0];
sx q[0];
rz(2.9184166) q[0];
rz(-3.0617833) q[1];
sx q[1];
rz(-2.1640919) q[1];
sx q[1];
rz(-0.83597437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30779926) q[0];
sx q[0];
rz(-0.5029486) q[0];
sx q[0];
rz(2.0849243) q[0];
x q[1];
rz(-0.67422949) q[2];
sx q[2];
rz(-2.5449356) q[2];
sx q[2];
rz(-0.7772738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8174396) q[1];
sx q[1];
rz(-2.0803703) q[1];
sx q[1];
rz(2.669599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23578819) q[3];
sx q[3];
rz(-0.97321586) q[3];
sx q[3];
rz(1.8862629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9072546) q[2];
sx q[2];
rz(-0.85504389) q[2];
sx q[2];
rz(-1.7503395) q[2];
rz(-1.2587345) q[3];
sx q[3];
rz(-1.3842868) q[3];
sx q[3];
rz(-2.7425227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626955) q[0];
sx q[0];
rz(-2.6628222) q[0];
sx q[0];
rz(0.71267772) q[0];
rz(-1.124294) q[1];
sx q[1];
rz(-1.5713888) q[1];
sx q[1];
rz(-1.0363151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881104) q[0];
sx q[0];
rz(-1.453791) q[0];
sx q[0];
rz(2.5036158) q[0];
rz(0.2377788) q[2];
sx q[2];
rz(-1.2847804) q[2];
sx q[2];
rz(-0.41508383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1829002) q[1];
sx q[1];
rz(-0.98366504) q[1];
sx q[1];
rz(-2.8927813) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18671496) q[3];
sx q[3];
rz(-2.2590504) q[3];
sx q[3];
rz(2.6783239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0816281) q[2];
sx q[2];
rz(-2.5689503) q[2];
sx q[2];
rz(-0.51679483) q[2];
rz(-1.9948657) q[3];
sx q[3];
rz(-2.1490993) q[3];
sx q[3];
rz(-3.0920046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.0741172) q[0];
sx q[0];
rz(-0.68454409) q[0];
sx q[0];
rz(-1.6790947) q[0];
rz(1.0985724) q[1];
sx q[1];
rz(-1.4417646) q[1];
sx q[1];
rz(-0.77902737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1027101) q[0];
sx q[0];
rz(-1.5691367) q[0];
sx q[0];
rz(3.0872905) q[0];
x q[1];
rz(-2.674496) q[2];
sx q[2];
rz(-1.0425241) q[2];
sx q[2];
rz(1.4548276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9866619) q[1];
sx q[1];
rz(-2.8720461) q[1];
sx q[1];
rz(2.4433072) q[1];
x q[2];
rz(-2.0059465) q[3];
sx q[3];
rz(-2.3823934) q[3];
sx q[3];
rz(-2.5569806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3940008) q[2];
sx q[2];
rz(-0.70351768) q[2];
sx q[2];
rz(-0.35959378) q[2];
rz(-2.8408585) q[3];
sx q[3];
rz(-0.81529236) q[3];
sx q[3];
rz(2.1451779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71036285) q[0];
sx q[0];
rz(-1.2082986) q[0];
sx q[0];
rz(2.8408458) q[0];
rz(2.6077479) q[1];
sx q[1];
rz(-2.0678935) q[1];
sx q[1];
rz(-0.11375443) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9886564) q[0];
sx q[0];
rz(-2.0242891) q[0];
sx q[0];
rz(-0.55983587) q[0];
x q[1];
rz(-3.0871307) q[2];
sx q[2];
rz(-0.60945933) q[2];
sx q[2];
rz(-0.0527339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2054322) q[1];
sx q[1];
rz(-1.1236842) q[1];
sx q[1];
rz(0.7128678) q[1];
x q[2];
rz(-1.9451009) q[3];
sx q[3];
rz(-1.7262207) q[3];
sx q[3];
rz(2.1799996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45812312) q[2];
sx q[2];
rz(-0.14294954) q[2];
sx q[2];
rz(0.45041034) q[2];
rz(-2.2266375) q[3];
sx q[3];
rz(-0.92792845) q[3];
sx q[3];
rz(-1.8042867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8704855) q[0];
sx q[0];
rz(-2.4021689) q[0];
sx q[0];
rz(-2.7258605) q[0];
rz(-0.43139002) q[1];
sx q[1];
rz(-2.5564671) q[1];
sx q[1];
rz(0.77692568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8857461) q[0];
sx q[0];
rz(-1.5579468) q[0];
sx q[0];
rz(1.5175411) q[0];
rz(2.326909) q[2];
sx q[2];
rz(-2.2534342) q[2];
sx q[2];
rz(-1.4614568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7207234) q[1];
sx q[1];
rz(-2.4178079) q[1];
sx q[1];
rz(-0.66026824) q[1];
rz(2.6747392) q[3];
sx q[3];
rz(-1.3131719) q[3];
sx q[3];
rz(0.063916884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1315769) q[2];
sx q[2];
rz(-2.3516529) q[2];
sx q[2];
rz(1.6291523) q[2];
rz(-0.43859628) q[3];
sx q[3];
rz(-2.3066543) q[3];
sx q[3];
rz(0.51630539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877391) q[0];
sx q[0];
rz(-2.6022311) q[0];
sx q[0];
rz(-1.3158276) q[0];
rz(-0.081789628) q[1];
sx q[1];
rz(-1.9606699) q[1];
sx q[1];
rz(1.2710424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91453856) q[0];
sx q[0];
rz(-2.7167121) q[0];
sx q[0];
rz(-1.6416613) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3009335) q[2];
sx q[2];
rz(-1.6191543) q[2];
sx q[2];
rz(-0.73634597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0473905) q[1];
sx q[1];
rz(-1.2823815) q[1];
sx q[1];
rz(-2.3408195) q[1];
x q[2];
rz(1.5958691) q[3];
sx q[3];
rz(-1.6524602) q[3];
sx q[3];
rz(1.9586399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3424993) q[2];
sx q[2];
rz(-1.4405595) q[2];
sx q[2];
rz(-1.052617) q[2];
rz(1.8140225) q[3];
sx q[3];
rz(-1.9460287) q[3];
sx q[3];
rz(-0.50935203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256182) q[0];
sx q[0];
rz(-1.7617891) q[0];
sx q[0];
rz(-1.2431385) q[0];
rz(-1.3432518) q[1];
sx q[1];
rz(-0.46011283) q[1];
sx q[1];
rz(-2.7759001) q[1];
rz(2.9467498) q[2];
sx q[2];
rz(-1.101782) q[2];
sx q[2];
rz(-0.89141104) q[2];
rz(1.5841019) q[3];
sx q[3];
rz(-1.9823488) q[3];
sx q[3];
rz(-2.3700598) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
