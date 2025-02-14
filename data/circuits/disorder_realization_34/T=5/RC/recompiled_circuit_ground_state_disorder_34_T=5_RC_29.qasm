OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(-2.5791383) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90949517) q[0];
sx q[0];
rz(-1.9350855) q[0];
sx q[0];
rz(-3.0845736) q[0];
rz(-pi) q[1];
rz(-0.15563528) q[2];
sx q[2];
rz(-2.0726207) q[2];
sx q[2];
rz(1.4841472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7592647) q[1];
sx q[1];
rz(-0.60324962) q[1];
sx q[1];
rz(-2.3690577) q[1];
rz(-pi) q[2];
rz(-0.39333087) q[3];
sx q[3];
rz(-1.1926023) q[3];
sx q[3];
rz(2.8761169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2354551) q[2];
sx q[2];
rz(-1.1780585) q[2];
sx q[2];
rz(-0.79929024) q[2];
rz(0.47131395) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039463194) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(1.348173) q[0];
rz(-1.7680291) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(1.9893533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593529) q[0];
sx q[0];
rz(-2.3639661) q[0];
sx q[0];
rz(-1.1471143) q[0];
x q[1];
rz(-0.67518465) q[2];
sx q[2];
rz(-2.084888) q[2];
sx q[2];
rz(-2.6509283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79346953) q[1];
sx q[1];
rz(-1.4033699) q[1];
sx q[1];
rz(-2.4504285) q[1];
x q[2];
rz(-0.51898414) q[3];
sx q[3];
rz(-2.0618084) q[3];
sx q[3];
rz(1.4903013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(0.93117923) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-2.0276766) q[3];
sx q[3];
rz(0.60025269) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057864144) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-0.51783836) q[0];
rz(0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(2.6944366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06181051) q[0];
sx q[0];
rz(-1.5070276) q[0];
sx q[0];
rz(0.26173862) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5340786) q[2];
sx q[2];
rz(-0.82582966) q[2];
sx q[2];
rz(-1.5811282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39147705) q[1];
sx q[1];
rz(-0.98408723) q[1];
sx q[1];
rz(0.16280414) q[1];
rz(-pi) q[2];
rz(2.6833862) q[3];
sx q[3];
rz(-2.5669006) q[3];
sx q[3];
rz(-1.1632869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7670224) q[2];
sx q[2];
rz(-0.76408237) q[2];
sx q[2];
rz(0.53263295) q[2];
rz(2.8940708) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(2.1029162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240876) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(-0.10661539) q[0];
rz(-0.34635776) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(1.1603629) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896987) q[0];
sx q[0];
rz(-1.3972613) q[0];
sx q[0];
rz(0.24994295) q[0];
x q[1];
rz(-0.71765064) q[2];
sx q[2];
rz(-0.96390488) q[2];
sx q[2];
rz(2.2362312) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44768279) q[1];
sx q[1];
rz(-0.47397787) q[1];
sx q[1];
rz(-0.86556566) q[1];
x q[2];
rz(-2.6909573) q[3];
sx q[3];
rz(-2.051578) q[3];
sx q[3];
rz(1.4522875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13566636) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(0.076210991) q[2];
rz(0.58586079) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(-1.7990254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0214486) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(-2.7914877) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(-1.9116481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6087225) q[0];
sx q[0];
rz(-1.2048831) q[0];
sx q[0];
rz(-2.8844112) q[0];
x q[1];
rz(-1.4054662) q[2];
sx q[2];
rz(-1.606719) q[2];
sx q[2];
rz(-1.7599811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38795162) q[1];
sx q[1];
rz(-1.3820096) q[1];
sx q[1];
rz(-1.2354047) q[1];
rz(-pi) q[2];
rz(1.2754945) q[3];
sx q[3];
rz(-0.60747889) q[3];
sx q[3];
rz(2.5584084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(2.7367075) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(0.83612061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(2.5591922) q[0];
rz(2.4549585) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.883421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67418881) q[0];
sx q[0];
rz(-1.6961401) q[0];
sx q[0];
rz(1.1320498) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35522009) q[2];
sx q[2];
rz(-2.0783011) q[2];
sx q[2];
rz(-1.2025646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2542779) q[1];
sx q[1];
rz(-1.5149024) q[1];
sx q[1];
rz(-3.0718551) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.889713) q[3];
sx q[3];
rz(-2.5237759) q[3];
sx q[3];
rz(-1.2417254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(-1.0180391) q[2];
rz(-2.8954519) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(1.6233981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(-2.5614118) q[0];
rz(-2.9980581) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-0.20763436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5136826) q[0];
sx q[0];
rz(-1.915689) q[0];
sx q[0];
rz(-0.78457997) q[0];
rz(-pi) q[1];
rz(2.0775547) q[2];
sx q[2];
rz(-2.4993976) q[2];
sx q[2];
rz(-2.3395777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5112111) q[1];
sx q[1];
rz(-2.2464754) q[1];
sx q[1];
rz(0.33501321) q[1];
rz(1.666388) q[3];
sx q[3];
rz(-1.2089653) q[3];
sx q[3];
rz(1.7887448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2988854) q[2];
sx q[2];
rz(-1.8536114) q[2];
sx q[2];
rz(0.60619727) q[2];
rz(-2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523478) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(1.0580753) q[0];
rz(-0.10969133) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.6995957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578975) q[0];
sx q[0];
rz(-0.51283118) q[0];
sx q[0];
rz(0.67016853) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8462366) q[2];
sx q[2];
rz(-1.0862964) q[2];
sx q[2];
rz(0.4415919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31046384) q[1];
sx q[1];
rz(-2.052124) q[1];
sx q[1];
rz(-0.83408611) q[1];
x q[2];
rz(0.85956802) q[3];
sx q[3];
rz(-0.35214927) q[3];
sx q[3];
rz(1.5052049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66120061) q[2];
sx q[2];
rz(-2.9471687) q[2];
sx q[2];
rz(1.2083758) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(-1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.1709568) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(-0.092967689) q[0];
rz(-1.2804821) q[1];
sx q[1];
rz(-0.72851506) q[1];
sx q[1];
rz(-0.13519898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8672392) q[0];
sx q[0];
rz(-1.6238191) q[0];
sx q[0];
rz(1.6163325) q[0];
rz(-pi) q[1];
rz(-0.5298631) q[2];
sx q[2];
rz(-2.8447897) q[2];
sx q[2];
rz(-0.57020818) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1982806) q[1];
sx q[1];
rz(-1.3168646) q[1];
sx q[1];
rz(-1.5928245) q[1];
x q[2];
rz(3.0673965) q[3];
sx q[3];
rz(-1.418494) q[3];
sx q[3];
rz(-2.0366813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70904237) q[2];
sx q[2];
rz(-1.9229527) q[2];
sx q[2];
rz(0.77152983) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563357) q[0];
sx q[0];
rz(-0.62194967) q[0];
sx q[0];
rz(1.1248032) q[0];
rz(-1.2987761) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-2.3769456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2128396) q[0];
sx q[0];
rz(-3.0055586) q[0];
sx q[0];
rz(-2.4555951) q[0];
rz(-pi) q[1];
rz(-2.6420399) q[2];
sx q[2];
rz(-2.1883114) q[2];
sx q[2];
rz(2.2411186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1435755) q[1];
sx q[1];
rz(-1.3680661) q[1];
sx q[1];
rz(0.12938093) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1922791) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(-1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(-0.20720227) q[2];
rz(2.1616705) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205263) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(1.4516713) q[2];
sx q[2];
rz(-1.4340192) q[2];
sx q[2];
rz(-1.7822969) q[2];
rz(-1.7332786) q[3];
sx q[3];
rz(-1.3357031) q[3];
sx q[3];
rz(1.6817844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
