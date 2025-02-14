OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0140822) q[0];
sx q[0];
rz(2.8988702) q[0];
sx q[0];
rz(10.299814) q[0];
rz(-1.8889282) q[1];
sx q[1];
rz(-1.6369605) q[1];
sx q[1];
rz(-0.59706444) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1738773) q[0];
sx q[0];
rz(-1.1861897) q[0];
sx q[0];
rz(-0.091716856) q[0];
rz(-pi) q[1];
rz(-0.70266415) q[2];
sx q[2];
rz(-1.4070373) q[2];
sx q[2];
rz(0.57715774) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7272612) q[1];
sx q[1];
rz(-2.4684069) q[1];
sx q[1];
rz(3.0607231) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4041599) q[3];
sx q[3];
rz(-1.5207411) q[3];
sx q[3];
rz(-1.6470294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.732932) q[2];
sx q[2];
rz(-0.084736846) q[2];
sx q[2];
rz(1.6049467) q[2];
rz(-1.599865) q[3];
sx q[3];
rz(-2.7432975) q[3];
sx q[3];
rz(-1.0140243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23564553) q[0];
sx q[0];
rz(-0.54536098) q[0];
sx q[0];
rz(-1.9769309) q[0];
rz(2.2230478) q[1];
sx q[1];
rz(-2.6530177) q[1];
sx q[1];
rz(-0.69893828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0889269) q[0];
sx q[0];
rz(-2.1667826) q[0];
sx q[0];
rz(1.2794897) q[0];
rz(-1.318803) q[2];
sx q[2];
rz(-2.4341035) q[2];
sx q[2];
rz(-1.8088264) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0731502) q[1];
sx q[1];
rz(-0.13769503) q[1];
sx q[1];
rz(-0.87611468) q[1];
x q[2];
rz(2.9925542) q[3];
sx q[3];
rz(-1.5939398) q[3];
sx q[3];
rz(1.4536957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32236448) q[2];
sx q[2];
rz(-1.4388204) q[2];
sx q[2];
rz(-1.9361852) q[2];
rz(2.9820005) q[3];
sx q[3];
rz(-1.7528088) q[3];
sx q[3];
rz(3.0157183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50652248) q[0];
sx q[0];
rz(-3.0634026) q[0];
sx q[0];
rz(3.0248094) q[0];
rz(2.0109239) q[1];
sx q[1];
rz(-0.10215452) q[1];
sx q[1];
rz(0.048361383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581586) q[0];
sx q[0];
rz(-1.6412705) q[0];
sx q[0];
rz(-1.6638046) q[0];
rz(-pi) q[1];
x q[1];
rz(1.618967) q[2];
sx q[2];
rz(-1.4496231) q[2];
sx q[2];
rz(0.19751422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54679388) q[1];
sx q[1];
rz(-0.76028667) q[1];
sx q[1];
rz(-2.0368963) q[1];
x q[2];
rz(-1.8675818) q[3];
sx q[3];
rz(-3.0078553) q[3];
sx q[3];
rz(-2.7406364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2241263) q[2];
sx q[2];
rz(-1.6240865) q[2];
sx q[2];
rz(-2.4833014) q[2];
rz(1.0607464) q[3];
sx q[3];
rz(-0.81183934) q[3];
sx q[3];
rz(-0.15061024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1310405) q[0];
sx q[0];
rz(-3.0578767) q[0];
sx q[0];
rz(0.87378275) q[0];
rz(-2.1281758) q[1];
sx q[1];
rz(-0.0031009379) q[1];
sx q[1];
rz(2.3168964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254941) q[0];
sx q[0];
rz(-1.6125729) q[0];
sx q[0];
rz(-3.1095563) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0295139) q[2];
sx q[2];
rz(-2.3397987) q[2];
sx q[2];
rz(0.0060334671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4383126) q[1];
sx q[1];
rz(-1.6973799) q[1];
sx q[1];
rz(1.8476608) q[1];
rz(-pi) q[2];
rz(-2.0788105) q[3];
sx q[3];
rz(-2.130878) q[3];
sx q[3];
rz(2.687272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1707477) q[2];
sx q[2];
rz(-2.8837995) q[2];
sx q[2];
rz(2.8802059) q[2];
rz(-1.1970674) q[3];
sx q[3];
rz(-1.7944929) q[3];
sx q[3];
rz(2.4813467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561387) q[0];
sx q[0];
rz(-0.11341299) q[0];
sx q[0];
rz(1.0979106) q[0];
rz(-0.55770355) q[1];
sx q[1];
rz(-0.028956078) q[1];
sx q[1];
rz(-1.2835693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52434873) q[0];
sx q[0];
rz(-1.6094406) q[0];
sx q[0];
rz(-1.6225553) q[0];
rz(-2.8800542) q[2];
sx q[2];
rz(-1.6439624) q[2];
sx q[2];
rz(-1.9436702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7456012) q[1];
sx q[1];
rz(-1.6874871) q[1];
sx q[1];
rz(0.50164739) q[1];
rz(2.5398272) q[3];
sx q[3];
rz(-2.2048809) q[3];
sx q[3];
rz(0.82179797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4881318) q[2];
sx q[2];
rz(-2.2296495) q[2];
sx q[2];
rz(-2.8690191) q[2];
rz(-1.9778947) q[3];
sx q[3];
rz(-1.3397763) q[3];
sx q[3];
rz(2.1986966) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96384752) q[0];
sx q[0];
rz(-0.073539428) q[0];
sx q[0];
rz(-2.9195926) q[0];
rz(1.47413) q[1];
sx q[1];
rz(-3.1167897) q[1];
sx q[1];
rz(-2.7892392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8604252) q[0];
sx q[0];
rz(-1.5727497) q[0];
sx q[0];
rz(-0.00020247799) q[0];
rz(-pi) q[1];
rz(0.84105305) q[2];
sx q[2];
rz(-1.111472) q[2];
sx q[2];
rz(-1.5265319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87599355) q[1];
sx q[1];
rz(-2.378646) q[1];
sx q[1];
rz(-2.250092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54103764) q[3];
sx q[3];
rz(-2.411708) q[3];
sx q[3];
rz(0.1950099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0148129) q[2];
sx q[2];
rz(-1.1819906) q[2];
sx q[2];
rz(-2.354055) q[2];
rz(3.0102503) q[3];
sx q[3];
rz(-2.1613439) q[3];
sx q[3];
rz(2.378715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545559) q[0];
sx q[0];
rz(-3.1058703) q[0];
sx q[0];
rz(-0.63294739) q[0];
rz(2.5542906) q[1];
sx q[1];
rz(-3.1226698) q[1];
sx q[1];
rz(2.6515554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5486262) q[0];
sx q[0];
rz(-1.6621309) q[0];
sx q[0];
rz(-1.7305602) q[0];
rz(-pi) q[1];
rz(2.7964726) q[2];
sx q[2];
rz(-1.39087) q[2];
sx q[2];
rz(-0.82054446) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1541802) q[1];
sx q[1];
rz(-1.9077717) q[1];
sx q[1];
rz(-0.67055871) q[1];
x q[2];
rz(2.1774749) q[3];
sx q[3];
rz(-2.5593966) q[3];
sx q[3];
rz(0.7403928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6106762) q[2];
sx q[2];
rz(-0.55675113) q[2];
sx q[2];
rz(-0.69182932) q[2];
rz(-2.1107215) q[3];
sx q[3];
rz(-0.94018799) q[3];
sx q[3];
rz(-0.33215365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0523025) q[0];
sx q[0];
rz(-0.060929935) q[0];
sx q[0];
rz(-2.6887509) q[0];
rz(2.7678658) q[1];
sx q[1];
rz(-3.0174297) q[1];
sx q[1];
rz(2.9492818) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34626353) q[0];
sx q[0];
rz(-1.6843616) q[0];
sx q[0];
rz(1.5862443) q[0];
rz(-pi) q[1];
rz(1.9779124) q[2];
sx q[2];
rz(-1.0958899) q[2];
sx q[2];
rz(-2.6559336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9050185) q[1];
sx q[1];
rz(-0.55779167) q[1];
sx q[1];
rz(1.2015427) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13914827) q[3];
sx q[3];
rz(-1.8915716) q[3];
sx q[3];
rz(-2.2389984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33334357) q[2];
sx q[2];
rz(-3.0875751) q[2];
sx q[2];
rz(2.5755908) q[2];
rz(-2.443215) q[3];
sx q[3];
rz(-1.799182) q[3];
sx q[3];
rz(-3.0799871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653387) q[0];
sx q[0];
rz(-2.7955671) q[0];
sx q[0];
rz(1.457343) q[0];
rz(2.1894042) q[1];
sx q[1];
rz(-3.1377276) q[1];
sx q[1];
rz(1.7239408) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3450556) q[0];
sx q[0];
rz(-0.36204189) q[0];
sx q[0];
rz(2.1276228) q[0];
rz(-pi) q[1];
rz(1.5308916) q[2];
sx q[2];
rz(-2.6171472) q[2];
sx q[2];
rz(0.1631069) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5385069) q[1];
sx q[1];
rz(-0.17854796) q[1];
sx q[1];
rz(0.18263765) q[1];
rz(-pi) q[2];
rz(1.438497) q[3];
sx q[3];
rz(-1.6198564) q[3];
sx q[3];
rz(2.8956312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62966043) q[2];
sx q[2];
rz(-0.28714445) q[2];
sx q[2];
rz(-0.49893898) q[2];
rz(2.387909) q[3];
sx q[3];
rz(-0.53871173) q[3];
sx q[3];
rz(0.57873571) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653487) q[0];
sx q[0];
rz(-0.25043273) q[0];
sx q[0];
rz(0.28755406) q[0];
rz(-2.4240671) q[1];
sx q[1];
rz(-3.0375807) q[1];
sx q[1];
rz(1.2984553) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037894) q[0];
sx q[0];
rz(-2.2713243) q[0];
sx q[0];
rz(-2.325969) q[0];
rz(-pi) q[1];
rz(-0.7602306) q[2];
sx q[2];
rz(-0.76334526) q[2];
sx q[2];
rz(1.476045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7804177) q[1];
sx q[1];
rz(-2.6589825) q[1];
sx q[1];
rz(0.9622099) q[1];
rz(-0.27661037) q[3];
sx q[3];
rz(-2.8146126) q[3];
sx q[3];
rz(0.17736926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8269877) q[2];
sx q[2];
rz(-1.8525476) q[2];
sx q[2];
rz(0.89024603) q[2];
rz(3.1173949) q[3];
sx q[3];
rz(-2.9399293) q[3];
sx q[3];
rz(0.81082398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.191962) q[0];
sx q[0];
rz(-2.5856615) q[0];
sx q[0];
rz(-0.91188201) q[0];
rz(0.98056071) q[1];
sx q[1];
rz(-2.801827) q[1];
sx q[1];
rz(-2.7580072) q[1];
rz(-1.8843731) q[2];
sx q[2];
rz(-1.6403392) q[2];
sx q[2];
rz(2.7210515) q[2];
rz(-2.2617658) q[3];
sx q[3];
rz(-0.36986989) q[3];
sx q[3];
rz(-1.5218745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
