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
rz(3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9026731) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(-0.33688776) q[0];
rz(-pi) q[1];
rz(0.84471976) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-0.5288045) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85841132) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(1.4285425) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48304708) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.7529863) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.7659448) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398867) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(0.4102104) q[0];
rz(-pi) q[1];
rz(-0.045250821) q[2];
sx q[2];
rz(-1.8088733) q[2];
sx q[2];
rz(1.0649875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34054204) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(1.3786475) q[1];
x q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.2352357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19757195) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(0.59528415) q[0];
rz(-pi) q[1];
rz(0.58358242) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(2.8180518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55360824) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(-1.8334465) q[1];
x q[2];
rz(-2.5123185) q[3];
sx q[3];
rz(-2.4089703) q[3];
sx q[3];
rz(-0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-3.0217357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170341) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(-1.6334565) q[0];
rz(-2.7718133) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-1.0880926) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82512292) q[1];
sx q[1];
rz(-2.9087062) q[1];
sx q[1];
rz(1.6474001) q[1];
rz(2.4265392) q[3];
sx q[3];
rz(-0.73521571) q[3];
sx q[3];
rz(2.8639567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.240085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86868984) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(2.9024283) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24057062) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(1.0742998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33927321) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(3.0250711) q[1];
rz(2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.6606768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(1.7646164) q[0];
rz(-pi) q[1];
rz(-0.00077498282) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(1.0926525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8112091) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.1793544) q[1];
rz(1.6310286) q[3];
sx q[3];
rz(-2.4619953) q[3];
sx q[3];
rz(-0.70590245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(-0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(3.0142473) q[0];
rz(0.84682805) q[2];
sx q[2];
rz(-0.70038751) q[2];
sx q[2];
rz(-0.91947739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8285671) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(-1.9810956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2277101) q[3];
sx q[3];
rz(-1.8808639) q[3];
sx q[3];
rz(0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(2.1929072) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(-0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.8519648) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-2.3873781) q[0];
sx q[0];
rz(-2.450599) q[0];
x q[1];
rz(2.3194359) q[2];
sx q[2];
rz(-0.90688721) q[2];
sx q[2];
rz(2.8234931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(-0.11282632) q[1];
rz(-pi) q[2];
rz(1.4054221) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(0.28265488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455505) q[0];
sx q[0];
rz(-2.1746785) q[0];
sx q[0];
rz(2.1349483) q[0];
x q[1];
rz(1.3547782) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(-0.87460364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.104056) q[1];
sx q[1];
rz(-1.9754585) q[1];
sx q[1];
rz(-1.7621103) q[1];
rz(-pi) q[2];
rz(-0.11741365) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(-2.5168602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012321) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(2.9451314) q[0];
rz(0.68306142) q[2];
sx q[2];
rz(-2.1944754) q[2];
sx q[2];
rz(1.1766441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1055923) q[1];
sx q[1];
rz(-2.3899374) q[1];
sx q[1];
rz(-0.73026258) q[1];
rz(-pi) q[2];
rz(-0.82716771) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(-3.0527715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
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
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.66424673) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(0.38469436) q[3];
sx q[3];
rz(-1.3833429) q[3];
sx q[3];
rz(-0.27237567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
