OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(-0.0012794415) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2389195) q[0];
sx q[0];
rz(-2.9905149) q[0];
sx q[0];
rz(-0.33688776) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(2.6127882) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85841132) q[1];
sx q[1];
rz(-1.6106669) q[1];
sx q[1];
rz(1.7130501) q[1];
x q[2];
rz(0.48304708) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.7529863) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028582024) q[0];
sx q[0];
rz(-0.41886371) q[0];
sx q[0];
rz(-0.21582614) q[0];
rz(-pi) q[1];
rz(-1.7550811) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(2.266303) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34054204) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(1.0248915) q[3];
sx q[3];
rz(-2.5149269) q[3];
sx q[3];
rz(2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.2352357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429758) q[0];
sx q[0];
rz(-1.9447864) q[0];
sx q[0];
rz(2.5234733) q[0];
rz(-pi) q[1];
rz(0.37249506) q[2];
sx q[2];
rz(-1.3349581) q[2];
sx q[2];
rz(-1.7861988) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8894549) q[1];
sx q[1];
rz(-0.368202) q[1];
sx q[1];
rz(-0.77106573) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5127605) q[3];
sx q[3];
rz(-1.1662081) q[3];
sx q[3];
rz(-1.7189327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45194295) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(-2.1769051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-1.0880926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4704628) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(-1.8030241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.8653387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.240085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2729028) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(-2.9024283) q[0];
x q[1];
rz(1.8334332) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(-2.5831985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8023194) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(3.0250711) q[1];
rz(-pi) q[2];
rz(-0.92091839) q[3];
sx q[3];
rz(-0.71483597) q[3];
sx q[3];
rz(-0.20540796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.038625) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.6606768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71125644) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(1.3769763) q[0];
x q[1];
rz(3.1408177) q[2];
sx q[2];
rz(-1.1256071) q[2];
sx q[2];
rz(-1.0926525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8112091) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.9622383) q[1];
rz(-pi) q[2];
rz(-3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(0.62852678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(-0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949188) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-0.17636756) q[0];
x q[1];
rz(-2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.2457459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8285671) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(1.9810956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38447325) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.2896279) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063110654) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(2.1102935) q[0];
x q[1];
rz(-2.3194359) q[2];
sx q[2];
rz(-0.90688721) q[2];
sx q[2];
rz(0.31809959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3370812) q[1];
sx q[1];
rz(-1.4599428) q[1];
sx q[1];
rz(-1.3827419) q[1];
rz(1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(-2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(-2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982272) q[0];
sx q[0];
rz(-0.8015612) q[0];
sx q[0];
rz(2.482224) q[0];
rz(-1.3547782) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(0.87460364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54284755) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(-0.41136841) q[1];
x q[2];
rz(-0.11741365) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(-2.5168602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1334575) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(2.125678) q[0];
x q[1];
rz(-2.4585312) q[2];
sx q[2];
rz(-0.94711727) q[2];
sx q[2];
rz(1.9649486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9224285) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(2.1283172) q[1];
rz(0.82716771) q[3];
sx q[3];
rz(-2.4359772) q[3];
sx q[3];
rz(0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-2.8005023) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-0.17102374) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-0.66424673) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(2.7568983) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
