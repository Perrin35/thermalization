OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(1.877797) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6440463) q[0];
sx q[0];
rz(-0.9263557) q[0];
sx q[0];
rz(-2.833616) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1793288) q[2];
sx q[2];
rz(-0.66676408) q[2];
sx q[2];
rz(-0.61362574) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25322029) q[1];
sx q[1];
rz(-1.651763) q[1];
sx q[1];
rz(-1.2451647) q[1];
x q[2];
rz(-2.0588082) q[3];
sx q[3];
rz(-2.6533439) q[3];
sx q[3];
rz(-1.4940718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15452142) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(-1.7729574) q[2];
rz(0.23855071) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(-3.0273048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40047613) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.1503295) q[0];
rz(0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(-1.5706583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99518665) q[0];
sx q[0];
rz(-1.0885462) q[0];
sx q[0];
rz(1.4508119) q[0];
x q[1];
rz(2.9628062) q[2];
sx q[2];
rz(-1.6331722) q[2];
sx q[2];
rz(1.5966001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1024061) q[1];
sx q[1];
rz(-0.79823179) q[1];
sx q[1];
rz(-0.95603966) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0474954) q[3];
sx q[3];
rz(-0.83752756) q[3];
sx q[3];
rz(2.1708716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34438434) q[2];
sx q[2];
rz(-2.4281561) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(2.3421085) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(-2.1605087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568273) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(0.10936603) q[0];
rz(-2.5378387) q[1];
sx q[1];
rz(-1.3394638) q[1];
sx q[1];
rz(-1.0341136) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84635705) q[0];
sx q[0];
rz(-0.24938008) q[0];
sx q[0];
rz(1.0037104) q[0];
x q[1];
rz(1.500152) q[2];
sx q[2];
rz(-1.6598668) q[2];
sx q[2];
rz(-1.8235109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2455403) q[1];
sx q[1];
rz(-2.8421092) q[1];
sx q[1];
rz(2.7761032) q[1];
x q[2];
rz(2.3699371) q[3];
sx q[3];
rz(-1.4691741) q[3];
sx q[3];
rz(1.9664498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5321396) q[2];
sx q[2];
rz(-2.1707462) q[2];
sx q[2];
rz(0.54507059) q[2];
rz(-2.3405781) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(-0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545749) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(2.6825478) q[0];
rz(2.0252939) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(-0.8265411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5640663) q[0];
sx q[0];
rz(-2.5788529) q[0];
sx q[0];
rz(2.2731645) q[0];
rz(-pi) q[1];
rz(1.6693833) q[2];
sx q[2];
rz(-1.4736097) q[2];
sx q[2];
rz(-2.5105385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53163995) q[1];
sx q[1];
rz(-2.5594257) q[1];
sx q[1];
rz(1.9332631) q[1];
x q[2];
rz(2.8119254) q[3];
sx q[3];
rz(-2.7894944) q[3];
sx q[3];
rz(2.6702777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35215846) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(0.74971548) q[2];
rz(2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4753251) q[0];
sx q[0];
rz(-2.1551977) q[0];
sx q[0];
rz(3.0295897) q[0];
rz(-2.9169967) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-0.78132838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3488779) q[0];
sx q[0];
rz(-1.9747866) q[0];
sx q[0];
rz(3.1125665) q[0];
rz(2.8728169) q[2];
sx q[2];
rz(-0.34160638) q[2];
sx q[2];
rz(0.20526055) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.427436) q[1];
sx q[1];
rz(-2.2793152) q[1];
sx q[1];
rz(-0.57517131) q[1];
rz(-2.9797793) q[3];
sx q[3];
rz(-2.1573632) q[3];
sx q[3];
rz(-2.2590421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3096699) q[2];
sx q[2];
rz(-2.2396542) q[2];
sx q[2];
rz(-2.208948) q[2];
rz(2.0617088) q[3];
sx q[3];
rz(-0.44411689) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90750736) q[0];
sx q[0];
rz(-2.0103173) q[0];
sx q[0];
rz(1.750741) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.876588) q[1];
sx q[1];
rz(-1.5914241) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54133139) q[0];
sx q[0];
rz(-2.2131767) q[0];
sx q[0];
rz(3.0976901) q[0];
rz(-3.0671547) q[2];
sx q[2];
rz(-0.46069579) q[2];
sx q[2];
rz(-1.4975394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6350123) q[1];
sx q[1];
rz(-1.6939347) q[1];
sx q[1];
rz(-1.413373) q[1];
x q[2];
rz(1.9523296) q[3];
sx q[3];
rz(-1.9100827) q[3];
sx q[3];
rz(-1.8041704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8288237) q[2];
sx q[2];
rz(-0.88461107) q[2];
sx q[2];
rz(-1.620232) q[2];
rz(-0.96794266) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(-2.2255285) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62436002) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(-2.970001) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(-1.4412057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65214163) q[0];
sx q[0];
rz(-1.2496557) q[0];
sx q[0];
rz(-1.7016861) q[0];
x q[1];
rz(-1.4588548) q[2];
sx q[2];
rz(-2.0840692) q[2];
sx q[2];
rz(-1.3789195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4413195) q[1];
sx q[1];
rz(-2.0236349) q[1];
sx q[1];
rz(-0.39703607) q[1];
rz(-pi) q[2];
rz(-3.0002241) q[3];
sx q[3];
rz(-1.0118359) q[3];
sx q[3];
rz(1.5593004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8695996) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(-0.54455152) q[2];
rz(0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(-0.47206363) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2917824) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(-2.612402) q[0];
rz(2.5657907) q[1];
sx q[1];
rz(-1.878783) q[1];
sx q[1];
rz(2.0226488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125583) q[0];
sx q[0];
rz(-0.055744113) q[0];
sx q[0];
rz(-2.7554465) q[0];
rz(-0.37092692) q[2];
sx q[2];
rz(-1.3998404) q[2];
sx q[2];
rz(1.8425187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0928438) q[1];
sx q[1];
rz(-2.1511937) q[1];
sx q[1];
rz(-1.012085) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20967926) q[3];
sx q[3];
rz(-2.0884313) q[3];
sx q[3];
rz(-0.28214327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.240856) q[2];
sx q[2];
rz(-0.46515981) q[2];
sx q[2];
rz(2.4294803) q[2];
rz(1.5322878) q[3];
sx q[3];
rz(-0.94523793) q[3];
sx q[3];
rz(1.3473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79494548) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(1.0711063) q[0];
rz(-1.935293) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(-1.4869022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3275571) q[0];
sx q[0];
rz(-2.2326755) q[0];
sx q[0];
rz(-1.8083841) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5096312) q[2];
sx q[2];
rz(-0.8425396) q[2];
sx q[2];
rz(-1.8441083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8739972) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(1.3895172) q[1];
rz(-pi) q[2];
rz(1.1219209) q[3];
sx q[3];
rz(-1.316183) q[3];
sx q[3];
rz(1.5790347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6335166) q[2];
sx q[2];
rz(-2.241892) q[2];
sx q[2];
rz(0.077795204) q[2];
rz(0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(-2.0859065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28854293) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(-0.37260923) q[0];
rz(0.44652069) q[1];
sx q[1];
rz(-1.6852854) q[1];
sx q[1];
rz(0.85652295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.646333) q[0];
sx q[0];
rz(-0.071852751) q[0];
sx q[0];
rz(2.4836579) q[0];
rz(2.2394453) q[2];
sx q[2];
rz(-2.9716316) q[2];
sx q[2];
rz(-0.68357498) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3831272) q[1];
sx q[1];
rz(-1.9680319) q[1];
sx q[1];
rz(-0.12786156) q[1];
x q[2];
rz(0.89337279) q[3];
sx q[3];
rz(-2.2489298) q[3];
sx q[3];
rz(-3.107389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(0.55553931) q[2];
rz(2.7555452) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1356708) q[0];
sx q[0];
rz(-1.4501403) q[0];
sx q[0];
rz(-1.8474664) q[0];
rz(-2.338943) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-0.31617185) q[2];
sx q[2];
rz(-2.7718622) q[2];
sx q[2];
rz(-2.3864701) q[2];
rz(0.84011806) q[3];
sx q[3];
rz(-1.6230604) q[3];
sx q[3];
rz(-0.75480672) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
