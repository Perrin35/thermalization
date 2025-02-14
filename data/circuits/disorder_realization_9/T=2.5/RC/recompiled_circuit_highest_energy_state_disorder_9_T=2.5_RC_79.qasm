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
rz(-5.4391556) q[1];
sx q[1];
rz(0.56517711) q[1];
sx q[1];
rz(11.204389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3484074) q[0];
sx q[0];
rz(-2.4688666) q[0];
sx q[0];
rz(2.5701218) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7704813) q[2];
sx q[2];
rz(-2.7058995) q[2];
sx q[2];
rz(-1.0267804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2842362) q[1];
sx q[1];
rz(-1.3487089) q[1];
sx q[1];
rz(2.7276993) q[1];
rz(-pi) q[2];
rz(0.99895044) q[3];
sx q[3];
rz(-2.5643454) q[3];
sx q[3];
rz(1.5791864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7094946) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(-2.438365) q[2];
rz(0.75913366) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(2.4659992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264939) q[0];
sx q[0];
rz(-2.4166985) q[0];
sx q[0];
rz(2.3469927) q[0];
rz(0.82851797) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(-0.4812831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3743418) q[0];
sx q[0];
rz(-1.6899979) q[0];
sx q[0];
rz(1.5425372) q[0];
x q[1];
rz(-2.0688624) q[2];
sx q[2];
rz(-0.8862555) q[2];
sx q[2];
rz(-1.7430151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6782234) q[1];
sx q[1];
rz(-1.3407523) q[1];
sx q[1];
rz(-0.078945625) q[1];
x q[2];
rz(-2.4614905) q[3];
sx q[3];
rz(-0.77723336) q[3];
sx q[3];
rz(-2.9297364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1543697) q[2];
sx q[2];
rz(-2.1798446) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(2.4842723) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.9067732) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364612) q[0];
sx q[0];
rz(-1.5345804) q[0];
sx q[0];
rz(2.8954647) q[0];
rz(-2.4283465) q[1];
sx q[1];
rz(-1.5245707) q[1];
sx q[1];
rz(2.2770142) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35219463) q[0];
sx q[0];
rz(-1.3653132) q[0];
sx q[0];
rz(1.8858389) q[0];
x q[1];
rz(0.53094338) q[2];
sx q[2];
rz(-2.1111672) q[2];
sx q[2];
rz(-0.8671538) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6990825) q[1];
sx q[1];
rz(-1.1103295) q[1];
sx q[1];
rz(0.090126474) q[1];
rz(-pi) q[2];
rz(-0.15154408) q[3];
sx q[3];
rz(-2.3484548) q[3];
sx q[3];
rz(2.0878937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41512179) q[2];
sx q[2];
rz(-2.7455726) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(-0.11780277) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(0.12731586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.4609818) q[0];
sx q[0];
rz(0.46517459) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(1.7101589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2003882) q[0];
sx q[0];
rz(-2.6058307) q[0];
sx q[0];
rz(0.8922116) q[0];
rz(-pi) q[1];
rz(-1.0776005) q[2];
sx q[2];
rz(-1.3921129) q[2];
sx q[2];
rz(0.57236949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5618973) q[1];
sx q[1];
rz(-2.2219147) q[1];
sx q[1];
rz(2.0148333) q[1];
x q[2];
rz(-1.9172098) q[3];
sx q[3];
rz(-2.2675121) q[3];
sx q[3];
rz(-0.55075607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(-1.6453936) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.8206785) q[3];
sx q[3];
rz(-2.3827609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.1372304) q[0];
sx q[0];
rz(-2.4431603) q[0];
sx q[0];
rz(1.6001562) q[0];
rz(1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.5024332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7961032) q[0];
sx q[0];
rz(-1.1123344) q[0];
sx q[0];
rz(2.084341) q[0];
rz(-pi) q[1];
rz(0.2404332) q[2];
sx q[2];
rz(-2.5747262) q[2];
sx q[2];
rz(2.4971003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1790601) q[1];
sx q[1];
rz(-1.1241364) q[1];
sx q[1];
rz(-2.0455749) q[1];
x q[2];
rz(-2.5859896) q[3];
sx q[3];
rz(-3.0137339) q[3];
sx q[3];
rz(3.1154261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0043103546) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(2.7610371) q[2];
rz(-2.6015094) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(1.9758457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49749097) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(0.38462001) q[0];
rz(3.0907471) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(1.6772038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.56154) q[0];
sx q[0];
rz(-1.1275575) q[0];
sx q[0];
rz(2.0361986) q[0];
x q[1];
rz(-2.5821677) q[2];
sx q[2];
rz(-0.96772268) q[2];
sx q[2];
rz(2.32351) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37300996) q[1];
sx q[1];
rz(-1.6554195) q[1];
sx q[1];
rz(0.090589295) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7767653) q[3];
sx q[3];
rz(-0.89650963) q[3];
sx q[3];
rz(-0.29335653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1579608) q[2];
sx q[2];
rz(-2.3775358) q[2];
sx q[2];
rz(0.91510406) q[2];
rz(0.30465952) q[3];
sx q[3];
rz(-2.4200078) q[3];
sx q[3];
rz(1.4753531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810998) q[0];
sx q[0];
rz(-2.2486794) q[0];
sx q[0];
rz(2.0258946) q[0];
rz(0.64104331) q[1];
sx q[1];
rz(-1.8843001) q[1];
sx q[1];
rz(1.3263652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64102252) q[0];
sx q[0];
rz(-1.4395073) q[0];
sx q[0];
rz(-0.48878756) q[0];
rz(-pi) q[1];
rz(0.37441476) q[2];
sx q[2];
rz(-1.4157989) q[2];
sx q[2];
rz(-2.5520476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5976482) q[1];
sx q[1];
rz(-0.94606384) q[1];
sx q[1];
rz(-2.2398857) q[1];
rz(0.85878813) q[3];
sx q[3];
rz(-2.1534216) q[3];
sx q[3];
rz(1.4020385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5850087) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(1.8741685) q[2];
rz(0.061138717) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1214445) q[0];
sx q[0];
rz(-2.3763438) q[0];
sx q[0];
rz(2.6159317) q[0];
rz(1.4875745) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(-2.9590327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199969) q[0];
sx q[0];
rz(-0.47843009) q[0];
sx q[0];
rz(2.6493376) q[0];
x q[1];
rz(-2.1514074) q[2];
sx q[2];
rz(-2.6117439) q[2];
sx q[2];
rz(-2.9030346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5044587) q[1];
sx q[1];
rz(-2.4877423) q[1];
sx q[1];
rz(-1.8549396) q[1];
rz(-2.8556999) q[3];
sx q[3];
rz(-2.0276311) q[3];
sx q[3];
rz(1.735812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(3.0328499) q[2];
rz(0.10793081) q[3];
sx q[3];
rz(-0.35522541) q[3];
sx q[3];
rz(1.7598553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050960798) q[0];
sx q[0];
rz(-1.2010295) q[0];
sx q[0];
rz(-2.4884124) q[0];
rz(1.8494891) q[1];
sx q[1];
rz(-0.2457681) q[1];
sx q[1];
rz(-0.076315708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298593) q[0];
sx q[0];
rz(-2.7329067) q[0];
sx q[0];
rz(0.9291533) q[0];
rz(-1.7387975) q[2];
sx q[2];
rz(-1.6036878) q[2];
sx q[2];
rz(-2.5891182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9477651) q[1];
sx q[1];
rz(-2.3725531) q[1];
sx q[1];
rz(0.2001708) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35450046) q[3];
sx q[3];
rz(-2.1700942) q[3];
sx q[3];
rz(-1.8219624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.896758) q[2];
sx q[2];
rz(-1.9568169) q[2];
sx q[2];
rz(-0.99311382) q[2];
rz(-1.3541597) q[3];
sx q[3];
rz(-2.6024151) q[3];
sx q[3];
rz(1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3058474) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(-2.848023) q[0];
rz(1.0319895) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(-0.31732496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.246884) q[0];
sx q[0];
rz(-2.5114759) q[0];
sx q[0];
rz(0.21960857) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34658771) q[2];
sx q[2];
rz(-2.5591597) q[2];
sx q[2];
rz(0.96894852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.369243) q[1];
sx q[1];
rz(-1.1398106) q[1];
sx q[1];
rz(-1.1915156) q[1];
rz(1.6280319) q[3];
sx q[3];
rz(-1.2737107) q[3];
sx q[3];
rz(3.0155268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62852922) q[2];
sx q[2];
rz(-2.8184012) q[2];
sx q[2];
rz(0.31488669) q[2];
rz(3.1079187) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(-1.4382188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.9411511) q[1];
sx q[1];
rz(-2.0574175) q[1];
sx q[1];
rz(2.5331694) q[1];
rz(3.0931531) q[2];
sx q[2];
rz(-1.5413399) q[2];
sx q[2];
rz(0.057609859) q[2];
rz(-0.55704883) q[3];
sx q[3];
rz(-1.7799041) q[3];
sx q[3];
rz(1.1064116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
