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
rz(3.0577793) q[0];
sx q[0];
rz(-0.35141355) q[0];
sx q[0];
rz(-2.0961528) q[0];
rz(-2.2975629) q[1];
sx q[1];
rz(-0.56517711) q[1];
sx q[1];
rz(-1.3619818) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3484074) q[0];
sx q[0];
rz(-2.4688666) q[0];
sx q[0];
rz(-0.5714709) q[0];
rz(-pi) q[1];
rz(-1.9988598) q[2];
sx q[2];
rz(-1.6546103) q[2];
sx q[2];
rz(0.36255896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9633544) q[1];
sx q[1];
rz(-0.46666086) q[1];
sx q[1];
rz(-2.6299824) q[1];
x q[2];
rz(1.0697738) q[3];
sx q[3];
rz(-1.27099) q[3];
sx q[3];
rz(-2.6553947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43209806) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(-2.438365) q[2];
rz(2.382459) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(0.67559344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21509875) q[0];
sx q[0];
rz(-2.4166985) q[0];
sx q[0];
rz(-2.3469927) q[0];
rz(0.82851797) q[1];
sx q[1];
rz(-1.0732032) q[1];
sx q[1];
rz(-2.6603096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3414087) q[0];
sx q[0];
rz(-1.5988549) q[0];
sx q[0];
rz(-3.0223439) q[0];
rz(0.52954604) q[2];
sx q[2];
rz(-2.3193619) q[2];
sx q[2];
rz(-2.1088496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6782234) q[1];
sx q[1];
rz(-1.8008403) q[1];
sx q[1];
rz(3.062647) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65297834) q[3];
sx q[3];
rz(-1.1140454) q[3];
sx q[3];
rz(-1.8819609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98722297) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(2.4842723) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(-1.2348194) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364612) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(-2.8954647) q[0];
rz(-0.71324619) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(2.2770142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789398) q[0];
sx q[0];
rz(-1.3653132) q[0];
sx q[0];
rz(1.2557538) q[0];
rz(-0.53094338) q[2];
sx q[2];
rz(-1.0304255) q[2];
sx q[2];
rz(2.2744389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6990825) q[1];
sx q[1];
rz(-1.1103295) q[1];
sx q[1];
rz(3.0514662) q[1];
rz(0.15154408) q[3];
sx q[3];
rz(-0.79313784) q[3];
sx q[3];
rz(-1.053699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41512179) q[2];
sx q[2];
rz(-2.7455726) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(-0.11780277) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(-3.0142768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.4609818) q[0];
sx q[0];
rz(0.46517459) q[0];
rz(0.90452114) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(-1.7101589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9412044) q[0];
sx q[0];
rz(-2.6058307) q[0];
sx q[0];
rz(2.2493811) q[0];
rz(-1.935237) q[2];
sx q[2];
rz(-0.52204692) q[2];
sx q[2];
rz(-2.4625157) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27181131) q[1];
sx q[1];
rz(-1.2220739) q[1];
sx q[1];
rz(-2.4407102) q[1];
rz(-pi) q[2];
rz(-1.2243829) q[3];
sx q[3];
rz(-0.87408057) q[3];
sx q[3];
rz(-0.55075607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(-1.6453936) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(2.3827609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0043623) q[0];
sx q[0];
rz(-2.4431603) q[0];
sx q[0];
rz(-1.5414365) q[0];
rz(-1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.6391594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3454895) q[0];
sx q[0];
rz(-2.0292583) q[0];
sx q[0];
rz(2.084341) q[0];
rz(-0.55372766) q[2];
sx q[2];
rz(-1.6990176) q[2];
sx q[2];
rz(2.011337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17325832) q[1];
sx q[1];
rz(-1.9957665) q[1];
sx q[1];
rz(2.6476184) q[1];
x q[2];
rz(2.5859896) q[3];
sx q[3];
rz(-0.1278588) q[3];
sx q[3];
rz(-0.026166548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0043103546) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(0.38055554) q[2];
rz(-2.6015094) q[3];
sx q[3];
rz(-1.8916847) q[3];
sx q[3];
rz(1.165747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6441017) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(-3.0907471) q[1];
sx q[1];
rz(-2.6652002) q[1];
sx q[1];
rz(-1.4643889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2028785) q[0];
sx q[0];
rz(-1.9882294) q[0];
sx q[0];
rz(2.6532618) q[0];
rz(-0.88846859) q[2];
sx q[2];
rz(-1.1184449) q[2];
sx q[2];
rz(-0.41146989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7685827) q[1];
sx q[1];
rz(-1.6554195) q[1];
sx q[1];
rz(0.090589295) q[1];
rz(-pi) q[2];
rz(-0.6847295) q[3];
sx q[3];
rz(-1.4103508) q[3];
sx q[3];
rz(-1.993865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9836318) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(-0.91510406) q[2];
rz(0.30465952) q[3];
sx q[3];
rz(-2.4200078) q[3];
sx q[3];
rz(-1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810998) q[0];
sx q[0];
rz(-0.89291328) q[0];
sx q[0];
rz(-2.0258946) q[0];
rz(-0.64104331) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(1.3263652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5005701) q[0];
sx q[0];
rz(-1.4395073) q[0];
sx q[0];
rz(-2.6528051) q[0];
x q[1];
rz(-0.4037598) q[2];
sx q[2];
rz(-0.40382622) q[2];
sx q[2];
rz(0.60690875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7515727) q[1];
sx q[1];
rz(-2.2604001) q[1];
sx q[1];
rz(2.431206) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2828045) q[3];
sx q[3];
rz(-2.1534216) q[3];
sx q[3];
rz(1.7395541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.556584) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(-1.8741685) q[2];
rz(3.0804539) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(0.88501969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0201482) q[0];
sx q[0];
rz(-0.76524884) q[0];
sx q[0];
rz(0.52566093) q[0];
rz(-1.4875745) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(-0.18255998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2935242) q[0];
sx q[0];
rz(-1.3514567) q[0];
sx q[0];
rz(2.7128986) q[0];
rz(0.31086773) q[2];
sx q[2];
rz(-1.1345049) q[2];
sx q[2];
rz(0.41149652) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2845738) q[1];
sx q[1];
rz(-2.1942687) q[1];
sx q[1];
rz(-0.2116043) q[1];
x q[2];
rz(-2.8556999) q[3];
sx q[3];
rz(-2.0276311) q[3];
sx q[3];
rz(-1.4057807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.93099) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(3.0328499) q[2];
rz(0.10793081) q[3];
sx q[3];
rz(-0.35522541) q[3];
sx q[3];
rz(-1.3817374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050960798) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(-0.65318024) q[0];
rz(1.2921035) q[1];
sx q[1];
rz(-0.2457681) q[1];
sx q[1];
rz(0.076315708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298593) q[0];
sx q[0];
rz(-0.40868592) q[0];
sx q[0];
rz(-0.9291533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7650899) q[2];
sx q[2];
rz(-0.17116088) q[2];
sx q[2];
rz(0.82679798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6725823) q[1];
sx q[1];
rz(-0.82090506) q[1];
sx q[1];
rz(1.3806848) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35450046) q[3];
sx q[3];
rz(-0.97149847) q[3];
sx q[3];
rz(-1.3196303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2448347) q[2];
sx q[2];
rz(-1.9568169) q[2];
sx q[2];
rz(2.1484788) q[2];
rz(-1.7874329) q[3];
sx q[3];
rz(-2.6024151) q[3];
sx q[3];
rz(1.4208992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3058474) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(2.848023) q[0];
rz(1.0319895) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(2.8242677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.246884) q[0];
sx q[0];
rz(-2.5114759) q[0];
sx q[0];
rz(-0.21960857) q[0];
rz(-2.5869682) q[2];
sx q[2];
rz(-1.3828424) q[2];
sx q[2];
rz(0.30890572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5342478) q[1];
sx q[1];
rz(-2.5754654) q[1];
sx q[1];
rz(-2.4636763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18472291) q[3];
sx q[3];
rz(-0.30238849) q[3];
sx q[3];
rz(-0.31935605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5130634) q[2];
sx q[2];
rz(-2.8184012) q[2];
sx q[2];
rz(0.31488669) q[2];
rz(-0.033673938) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2122129) q[0];
sx q[0];
rz(-1.0502945) q[0];
sx q[0];
rz(2.1631277) q[0];
rz(-0.20044151) q[1];
sx q[1];
rz(-1.0841752) q[1];
sx q[1];
rz(-0.60842327) q[1];
rz(0.048439518) q[2];
sx q[2];
rz(-1.6002527) q[2];
sx q[2];
rz(-3.0839828) q[2];
rz(-2.5845438) q[3];
sx q[3];
rz(-1.3616886) q[3];
sx q[3];
rz(-2.035181) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
