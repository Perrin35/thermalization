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
rz(1.053831) q[0];
sx q[0];
rz(3.6917917) q[0];
sx q[0];
rz(10.79296) q[0];
rz(2.6693681) q[1];
sx q[1];
rz(-0.11657403) q[1];
sx q[1];
rz(0.20224686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41101) q[0];
sx q[0];
rz(-2.8484593) q[0];
sx q[0];
rz(2.3877451) q[0];
x q[1];
rz(-2.6361385) q[2];
sx q[2];
rz(-1.1040282) q[2];
sx q[2];
rz(2.9708178) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7649012) q[1];
sx q[1];
rz(-2.4623532) q[1];
sx q[1];
rz(2.0308073) q[1];
x q[2];
rz(1.1047045) q[3];
sx q[3];
rz(-1.6045609) q[3];
sx q[3];
rz(-2.1347783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5998485) q[2];
sx q[2];
rz(-2.2969963) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(0.26722515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32352725) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(0.7793119) q[0];
rz(1.3620954) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(2.0077226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529566) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(3.033328) q[0];
rz(-pi) q[1];
rz(-1.629584) q[2];
sx q[2];
rz(-1.4530621) q[2];
sx q[2];
rz(1.1716154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.081959978) q[1];
sx q[1];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.3038941) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0615223) q[3];
sx q[3];
rz(-1.5716553) q[3];
sx q[3];
rz(0.96405503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4506932) q[2];
sx q[2];
rz(-1.0319812) q[2];
sx q[2];
rz(0.1965941) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.5195547) q[3];
sx q[3];
rz(1.2381747) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152385) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(-0.73177904) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-2.7751353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55933773) q[0];
sx q[0];
rz(-0.32024239) q[0];
sx q[0];
rz(-2.5473464) q[0];
x q[1];
rz(-1.7973034) q[2];
sx q[2];
rz(-1.8069805) q[2];
sx q[2];
rz(0.073848595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4225821) q[1];
sx q[1];
rz(-1.9937464) q[1];
sx q[1];
rz(-0.99653901) q[1];
rz(-pi) q[2];
rz(2.2711804) q[3];
sx q[3];
rz(-1.0910209) q[3];
sx q[3];
rz(1.9241672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(2.2491573) q[2];
rz(2.6871032) q[3];
sx q[3];
rz(-0.82388866) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3733805) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(-0.40661231) q[0];
rz(-0.82798249) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(-0.83555317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53518772) q[0];
sx q[0];
rz(-1.7514896) q[0];
sx q[0];
rz(-0.040099696) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2242975) q[2];
sx q[2];
rz(-1.3801127) q[2];
sx q[2];
rz(0.71587038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20488508) q[1];
sx q[1];
rz(-0.46006535) q[1];
sx q[1];
rz(-2.8060444) q[1];
rz(-pi) q[2];
rz(0.097483272) q[3];
sx q[3];
rz(-1.1872224) q[3];
sx q[3];
rz(0.95065476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(1.3678) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(-0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-2.9827523) q[0];
rz(-1.5171492) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(-1.3295757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8027199) q[0];
sx q[0];
rz(-3.1204528) q[0];
sx q[0];
rz(2.6129524) q[0];
rz(-pi) q[1];
rz(-2.9022568) q[2];
sx q[2];
rz(-2.1194601) q[2];
sx q[2];
rz(1.2557097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6318887) q[1];
sx q[1];
rz(-1.77138) q[1];
sx q[1];
rz(1.9575809) q[1];
x q[2];
rz(-2.9333889) q[3];
sx q[3];
rz(-2.1284687) q[3];
sx q[3];
rz(-2.4488136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(-2.5895183) q[2];
rz(-2.3479346) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119222) q[0];
sx q[0];
rz(-0.86450082) q[0];
sx q[0];
rz(0.70575869) q[0];
rz(0.13748473) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(-0.71298832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891083) q[0];
sx q[0];
rz(-0.45408598) q[0];
sx q[0];
rz(-0.54291351) q[0];
rz(-2.627264) q[2];
sx q[2];
rz(-1.9514958) q[2];
sx q[2];
rz(-1.9388388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9116152) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(-1.9703034) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8637795) q[3];
sx q[3];
rz(-1.4594541) q[3];
sx q[3];
rz(-0.35929832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9299499) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(2.8575274) q[2];
rz(3.0197213) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26509869) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(-2.4342243) q[0];
rz(2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.7792938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.8168983) q[0];
sx q[0];
rz(-1.2854693) q[0];
x q[1];
rz(0.84378924) q[2];
sx q[2];
rz(-1.2449045) q[2];
sx q[2];
rz(1.0032723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55864298) q[1];
sx q[1];
rz(-2.042965) q[1];
sx q[1];
rz(-2.1583945) q[1];
rz(0.23162095) q[3];
sx q[3];
rz(-0.4670139) q[3];
sx q[3];
rz(-1.8927595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1661487) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-2.2834957) q[3];
sx q[3];
rz(-0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910832) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(-1.9367223) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(-2.9606294) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807225) q[0];
sx q[0];
rz(-1.1747169) q[0];
sx q[0];
rz(-1.164308) q[0];
rz(-1.7082735) q[2];
sx q[2];
rz(-1.2243234) q[2];
sx q[2];
rz(2.8203143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6540203) q[1];
sx q[1];
rz(-1.0165689) q[1];
sx q[1];
rz(0.2481064) q[1];
x q[2];
rz(-0.11060235) q[3];
sx q[3];
rz(-1.4975274) q[3];
sx q[3];
rz(-1.4721118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5391431) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(-0.55475956) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(0.35840148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48034126) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(3.0338147) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(1.617618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.658889) q[0];
sx q[0];
rz(-1.802717) q[0];
sx q[0];
rz(0.13302444) q[0];
rz(-pi) q[1];
rz(2.4977106) q[2];
sx q[2];
rz(-1.1538299) q[2];
sx q[2];
rz(1.8089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0131993) q[1];
sx q[1];
rz(-0.46401641) q[1];
sx q[1];
rz(1.6676519) q[1];
rz(-pi) q[2];
rz(0.49160853) q[3];
sx q[3];
rz(-1.5244686) q[3];
sx q[3];
rz(0.87418782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9109351) q[2];
sx q[2];
rz(-2.6378938) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(2.3835278) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199353) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(0.94104952) q[0];
rz(2.7811116) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(-1.2233268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62241232) q[0];
sx q[0];
rz(-2.522058) q[0];
sx q[0];
rz(-0.34468083) q[0];
rz(0.98127301) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(-1.5379932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66050038) q[1];
sx q[1];
rz(-1.6410488) q[1];
sx q[1];
rz(-0.41892799) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1072708) q[3];
sx q[3];
rz(-1.7400898) q[3];
sx q[3];
rz(-0.062549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2337522) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(0.20445538) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37888708) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(-1.0574404) q[2];
sx q[2];
rz(-2.3427137) q[2];
sx q[2];
rz(-2.8165934) q[2];
rz(2.8818535) q[3];
sx q[3];
rz(-2.9208214) q[3];
sx q[3];
rz(-0.95820273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
