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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(0.52752703) q[1];
sx q[1];
rz(-2.1105284) q[1];
sx q[1];
rz(0.48479015) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82370347) q[0];
sx q[0];
rz(-1.4394338) q[0];
sx q[0];
rz(-2.0844368) q[0];
x q[1];
rz(-0.77930696) q[2];
sx q[2];
rz(-1.8911084) q[2];
sx q[2];
rz(2.8576098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5875549) q[1];
sx q[1];
rz(-1.629519) q[1];
sx q[1];
rz(-2.1802933) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6438585) q[3];
sx q[3];
rz(-1.2803852) q[3];
sx q[3];
rz(2.6892457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.506044) q[2];
sx q[2];
rz(-0.94749331) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(-0.55073589) q[3];
sx q[3];
rz(-2.4516055) q[3];
sx q[3];
rz(-2.0042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15744844) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(2.3781811) q[0];
rz(-2.5570671) q[1];
sx q[1];
rz(-1.7598563) q[1];
sx q[1];
rz(2.1302628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69788591) q[0];
sx q[0];
rz(-2.9666565) q[0];
sx q[0];
rz(-0.849443) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23638921) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(-0.37416247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71537575) q[1];
sx q[1];
rz(-1.9131887) q[1];
sx q[1];
rz(-2.9222548) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743558) q[3];
sx q[3];
rz(-1.0938489) q[3];
sx q[3];
rz(-1.2881035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4804907) q[2];
sx q[2];
rz(-1.2031518) q[2];
sx q[2];
rz(2.1413546) q[2];
rz(-0.58146042) q[3];
sx q[3];
rz(-0.73271078) q[3];
sx q[3];
rz(1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30705273) q[0];
sx q[0];
rz(-2.5426799) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(0.3071951) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(-2.6934326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7355928) q[0];
sx q[0];
rz(-1.7032529) q[0];
sx q[0];
rz(2.622582) q[0];
x q[1];
rz(1.0501625) q[2];
sx q[2];
rz(-1.3853645) q[2];
sx q[2];
rz(2.7523169) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.243147) q[1];
sx q[1];
rz(-2.2013651) q[1];
sx q[1];
rz(3.10792) q[1];
rz(2.7072315) q[3];
sx q[3];
rz(-1.5944527) q[3];
sx q[3];
rz(2.0659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3192516) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(-0.59784293) q[2];
rz(2.6939825) q[3];
sx q[3];
rz(-1.8748583) q[3];
sx q[3];
rz(0.058535695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41281259) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(-0.70984167) q[0];
rz(2.1624883) q[1];
sx q[1];
rz(-0.85936463) q[1];
sx q[1];
rz(-1.8202579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036572) q[0];
sx q[0];
rz(-2.4383579) q[0];
sx q[0];
rz(-0.08589311) q[0];
x q[1];
rz(-2.3611907) q[2];
sx q[2];
rz(-1.4596887) q[2];
sx q[2];
rz(1.7826796) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7497341) q[1];
sx q[1];
rz(-0.94288153) q[1];
sx q[1];
rz(-2.8086752) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0274326) q[3];
sx q[3];
rz(-1.8413944) q[3];
sx q[3];
rz(0.71934127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1506302) q[2];
sx q[2];
rz(-1.4272828) q[2];
sx q[2];
rz(0.70685753) q[2];
rz(-1.3458716) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(-0.70485392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83955806) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(1.5020405) q[0];
rz(1.543965) q[1];
sx q[1];
rz(-2.2815956) q[1];
sx q[1];
rz(1.4942716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2199545) q[0];
sx q[0];
rz(-1.0003547) q[0];
sx q[0];
rz(3.1396237) q[0];
rz(-pi) q[1];
rz(0.59539184) q[2];
sx q[2];
rz(-1.5768345) q[2];
sx q[2];
rz(1.572682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9372796) q[1];
sx q[1];
rz(-1.2741852) q[1];
sx q[1];
rz(1.2526172) q[1];
rz(-0.022541209) q[3];
sx q[3];
rz(-2.1858741) q[3];
sx q[3];
rz(-0.32451567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8813701) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(2.5831045) q[2];
rz(-0.50885606) q[3];
sx q[3];
rz(-1.6130092) q[3];
sx q[3];
rz(1.8744899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5452071) q[0];
sx q[0];
rz(-1.8983497) q[0];
sx q[0];
rz(0.16045706) q[0];
rz(-2.4682553) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(1.1268667) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50235275) q[0];
sx q[0];
rz(-1.9316088) q[0];
sx q[0];
rz(-0.32397425) q[0];
rz(-pi) q[1];
rz(-2.3773002) q[2];
sx q[2];
rz(-1.8074161) q[2];
sx q[2];
rz(-2.3911084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4194132) q[1];
sx q[1];
rz(-0.99687885) q[1];
sx q[1];
rz(-2.9001321) q[1];
rz(-2.1127719) q[3];
sx q[3];
rz(-1.3888956) q[3];
sx q[3];
rz(1.0396837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1274073) q[2];
sx q[2];
rz(-1.5994453) q[2];
sx q[2];
rz(-0.7737774) q[2];
rz(2.1524147) q[3];
sx q[3];
rz(-2.7343605) q[3];
sx q[3];
rz(-1.0729084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147375) q[0];
sx q[0];
rz(-1.6586774) q[0];
sx q[0];
rz(-2.3684655) q[0];
rz(1.585656) q[1];
sx q[1];
rz(-1.8525886) q[1];
sx q[1];
rz(-1.1770491) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0814514) q[0];
sx q[0];
rz(-0.8340652) q[0];
sx q[0];
rz(3.0990776) q[0];
rz(-0.3118722) q[2];
sx q[2];
rz(-1.0200899) q[2];
sx q[2];
rz(0.68974173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29876626) q[1];
sx q[1];
rz(-2.2027059) q[1];
sx q[1];
rz(-0.34159391) q[1];
rz(-pi) q[2];
rz(-1.4312533) q[3];
sx q[3];
rz(-1.3948895) q[3];
sx q[3];
rz(1.6493451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3561463) q[2];
sx q[2];
rz(-0.36474228) q[2];
sx q[2];
rz(2.6981603) q[2];
rz(-2.7555079) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(0.20080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3767553) q[0];
sx q[0];
rz(-1.5489464) q[0];
sx q[0];
rz(-3.1078597) q[0];
rz(-2.4434166) q[1];
sx q[1];
rz(-0.59711421) q[1];
sx q[1];
rz(2.4765292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2584943) q[0];
sx q[0];
rz(-1.3083212) q[0];
sx q[0];
rz(-1.1061944) q[0];
rz(-pi) q[1];
rz(0.4735293) q[2];
sx q[2];
rz(-0.90190926) q[2];
sx q[2];
rz(-2.7471126) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57369697) q[1];
sx q[1];
rz(-2.0665281) q[1];
sx q[1];
rz(-2.2894409) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0785776) q[3];
sx q[3];
rz(-0.32013963) q[3];
sx q[3];
rz(-0.053016114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4195121) q[2];
sx q[2];
rz(-2.2245912) q[2];
sx q[2];
rz(1.3631932) q[2];
rz(-0.52115399) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1727961) q[0];
sx q[0];
rz(-1.3461312) q[0];
sx q[0];
rz(2.3222493) q[0];
rz(-2.4567545) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(-3.0192764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8259895) q[0];
sx q[0];
rz(-1.9902744) q[0];
sx q[0];
rz(-0.88780888) q[0];
x q[1];
rz(-2.165602) q[2];
sx q[2];
rz(-1.8907818) q[2];
sx q[2];
rz(0.94727635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6284991) q[1];
sx q[1];
rz(-0.6969531) q[1];
sx q[1];
rz(-2.6234187) q[1];
rz(-pi) q[2];
rz(1.1055533) q[3];
sx q[3];
rz(-2.1085582) q[3];
sx q[3];
rz(0.28013438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.704432) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(-0.86755794) q[2];
rz(-1.3982754) q[3];
sx q[3];
rz(-0.38845348) q[3];
sx q[3];
rz(0.5790264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.687998) q[0];
sx q[0];
rz(-1.9168251) q[0];
sx q[0];
rz(-2.0523742) q[0];
rz(3.094063) q[1];
sx q[1];
rz(-1.0462953) q[1];
sx q[1];
rz(-0.68738031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7318667) q[0];
sx q[0];
rz(-2.0610524) q[0];
sx q[0];
rz(-0.027564647) q[0];
rz(-pi) q[1];
rz(2.0014761) q[2];
sx q[2];
rz(-1.8269369) q[2];
sx q[2];
rz(1.3821951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7043651) q[1];
sx q[1];
rz(-0.95672551) q[1];
sx q[1];
rz(-1.006587) q[1];
x q[2];
rz(3.0484588) q[3];
sx q[3];
rz(-0.93619213) q[3];
sx q[3];
rz(0.30090162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99531949) q[2];
sx q[2];
rz(-2.5226888) q[2];
sx q[2];
rz(-1.4825561) q[2];
rz(-2.4847374) q[3];
sx q[3];
rz(-1.1494136) q[3];
sx q[3];
rz(-2.759554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8889846) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(-2.5629015) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(-0.22004057) q[2];
sx q[2];
rz(-0.36424988) q[2];
sx q[2];
rz(3.0314432) q[2];
rz(-0.5169576) q[3];
sx q[3];
rz(-1.4993954) q[3];
sx q[3];
rz(-2.7837337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
