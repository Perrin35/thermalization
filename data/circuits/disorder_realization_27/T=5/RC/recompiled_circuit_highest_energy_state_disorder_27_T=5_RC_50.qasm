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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(-2.9297096) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402474) q[0];
sx q[0];
rz(-1.5703526) q[0];
sx q[0];
rz(1.5612649) q[0];
x q[1];
rz(1.9340408) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(-2.5763047) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4204882) q[1];
sx q[1];
rz(-1.8454362) q[1];
sx q[1];
rz(0.65914776) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7805232) q[3];
sx q[3];
rz(-0.98376432) q[3];
sx q[3];
rz(-0.39654708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2892896) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(-2.1763109) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(0.088951237) q[0];
rz(1.7695919) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7331284) q[0];
sx q[0];
rz(-0.95703546) q[0];
sx q[0];
rz(-1.2888786) q[0];
rz(-pi) q[1];
rz(0.78271336) q[2];
sx q[2];
rz(-1.9316976) q[2];
sx q[2];
rz(-0.95907839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1762878) q[1];
sx q[1];
rz(-1.7200279) q[1];
sx q[1];
rz(2.01841) q[1];
rz(-pi) q[2];
rz(-2.6022909) q[3];
sx q[3];
rz(-1.4852583) q[3];
sx q[3];
rz(-2.364769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8770807) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(1.6748927) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(-0.27798852) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-0.40036449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63567296) q[0];
sx q[0];
rz(-2.189496) q[0];
sx q[0];
rz(1.0158422) q[0];
rz(-pi) q[1];
rz(0.30817356) q[2];
sx q[2];
rz(-0.71012596) q[2];
sx q[2];
rz(-1.9911204) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0644933) q[1];
sx q[1];
rz(-1.1005741) q[1];
sx q[1];
rz(-1.3920067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2648029) q[3];
sx q[3];
rz(-1.8974432) q[3];
sx q[3];
rz(1.5375002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(2.686783) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95551816) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(-3.1410826) q[0];
rz(0.60091248) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(-2.9972163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3568319) q[0];
sx q[0];
rz(-2.5481173) q[0];
sx q[0];
rz(-1.9463825) q[0];
rz(-pi) q[1];
rz(-0.52835502) q[2];
sx q[2];
rz(-2.8502185) q[2];
sx q[2];
rz(2.0896512) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3444654) q[1];
sx q[1];
rz(-0.68736156) q[1];
sx q[1];
rz(-2.7784154) q[1];
x q[2];
rz(-2.5944806) q[3];
sx q[3];
rz(-1.662385) q[3];
sx q[3];
rz(-1.6325337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(-1.9566253) q[0];
rz(-2.9153337) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(-1.0386946) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3861235) q[0];
sx q[0];
rz(-0.9842819) q[0];
sx q[0];
rz(1.6688132) q[0];
x q[1];
rz(-2.4034385) q[2];
sx q[2];
rz(-2.0655491) q[2];
sx q[2];
rz(-0.02502266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2782368) q[1];
sx q[1];
rz(-2.0788631) q[1];
sx q[1];
rz(0.79973508) q[1];
rz(1.1306954) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(-1.1820933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(-2.8254438) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(2.5441817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53086583) q[0];
sx q[0];
rz(-1.7590176) q[0];
sx q[0];
rz(1.782531) q[0];
x q[1];
rz(-2.8812203) q[2];
sx q[2];
rz(-0.7625167) q[2];
sx q[2];
rz(-2.0338361) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0874071) q[1];
sx q[1];
rz(-0.30255908) q[1];
sx q[1];
rz(2.5564479) q[1];
x q[2];
rz(3.0496063) q[3];
sx q[3];
rz(-1.2013271) q[3];
sx q[3];
rz(-1.8487768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(1.0763947) q[2];
rz(-1.9988029) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-3.1031188) q[0];
rz(-3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(2.9077392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215639) q[0];
sx q[0];
rz(-1.3498684) q[0];
sx q[0];
rz(1.7608282) q[0];
x q[1];
rz(1.2042768) q[2];
sx q[2];
rz(-1.8631753) q[2];
sx q[2];
rz(-0.93380837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5906787) q[1];
sx q[1];
rz(-1.6603025) q[1];
sx q[1];
rz(0.64915912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6285628) q[3];
sx q[3];
rz(-1.4411297) q[3];
sx q[3];
rz(2.4313789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-2.5433507) q[2];
rz(3.0277142) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(0.2463499) q[0];
rz(1.2975289) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(0.49682239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(2.4606649) q[0];
rz(-pi) q[1];
rz(-0.70513983) q[2];
sx q[2];
rz(-1.6821096) q[2];
sx q[2];
rz(1.706858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(-1.103765) q[1];
rz(-pi) q[2];
rz(-1.0786177) q[3];
sx q[3];
rz(-1.7075305) q[3];
sx q[3];
rz(1.1059974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.9971087) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(-2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(1.9937531) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(0.60417169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2691374) q[0];
sx q[0];
rz(-0.34664422) q[0];
sx q[0];
rz(0.64328648) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8762875) q[2];
sx q[2];
rz(-2.5992706) q[2];
sx q[2];
rz(3.0408183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3499602) q[1];
sx q[1];
rz(-0.59795982) q[1];
sx q[1];
rz(2.7954742) q[1];
rz(2.484453) q[3];
sx q[3];
rz(-0.98987386) q[3];
sx q[3];
rz(2.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26312795) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(-1.912502) q[3];
sx q[3];
rz(-1.2058328) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1935254) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(-2.926631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8037655) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(2.4830677) q[0];
x q[1];
rz(-2.6097492) q[2];
sx q[2];
rz(-1.1924679) q[2];
sx q[2];
rz(1.7500306) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4196288) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(2.2339905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39542012) q[3];
sx q[3];
rz(-0.59874615) q[3];
sx q[3];
rz(1.0637763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(-2.1827533) q[2];
rz(1.1622608) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6017629) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(-2.4330347) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(1.2461927) q[2];
sx q[2];
rz(-2.1094252) q[2];
sx q[2];
rz(-0.88627041) q[2];
rz(-2.419653) q[3];
sx q[3];
rz(-1.4293213) q[3];
sx q[3];
rz(-2.8534129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
