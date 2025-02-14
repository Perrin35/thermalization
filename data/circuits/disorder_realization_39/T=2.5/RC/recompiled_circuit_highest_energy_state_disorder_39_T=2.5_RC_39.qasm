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
rz(1.0562309) q[0];
sx q[0];
rz(-2.8662999) q[0];
sx q[0];
rz(-0.054962602) q[0];
rz(0.7420525) q[1];
sx q[1];
rz(3.2522178) q[1];
sx q[1];
rz(10.179959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7464712) q[0];
sx q[0];
rz(-2.3262505) q[0];
sx q[0];
rz(-1.59207) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58671592) q[2];
sx q[2];
rz(-2.3745076) q[2];
sx q[2];
rz(1.6471582) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29355907) q[1];
sx q[1];
rz(-1.5747442) q[1];
sx q[1];
rz(1.1396386) q[1];
rz(-pi) q[2];
rz(-2.8922454) q[3];
sx q[3];
rz(-0.60725242) q[3];
sx q[3];
rz(-0.23199122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5893958) q[2];
sx q[2];
rz(-1.9856717) q[2];
sx q[2];
rz(-1.7414306) q[2];
rz(-0.33341148) q[3];
sx q[3];
rz(-0.64781487) q[3];
sx q[3];
rz(-2.6441372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9452776) q[0];
sx q[0];
rz(-0.59756398) q[0];
sx q[0];
rz(3.0481098) q[0];
rz(-0.64158332) q[1];
sx q[1];
rz(-0.43951324) q[1];
sx q[1];
rz(0.19215597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963325) q[0];
sx q[0];
rz(-2.30034) q[0];
sx q[0];
rz(-1.6044046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21382217) q[2];
sx q[2];
rz(-1.3665843) q[2];
sx q[2];
rz(1.7776877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9201492) q[1];
sx q[1];
rz(-2.0196947) q[1];
sx q[1];
rz(-0.68660983) q[1];
rz(-pi) q[2];
x q[2];
rz(0.014161807) q[3];
sx q[3];
rz(-2.8175655) q[3];
sx q[3];
rz(-0.06122804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41500348) q[2];
sx q[2];
rz(-0.64565349) q[2];
sx q[2];
rz(-1.0059221) q[2];
rz(3.1065324) q[3];
sx q[3];
rz(-0.56060767) q[3];
sx q[3];
rz(0.75986552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1890892) q[0];
sx q[0];
rz(-1.3409505) q[0];
sx q[0];
rz(2.9553318) q[0];
rz(-0.26178905) q[1];
sx q[1];
rz(-1.9354542) q[1];
sx q[1];
rz(0.50938767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88424524) q[0];
sx q[0];
rz(-1.5871947) q[0];
sx q[0];
rz(-0.00014033982) q[0];
rz(-0.89318067) q[2];
sx q[2];
rz(-1.2959157) q[2];
sx q[2];
rz(-2.6468697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2434329) q[1];
sx q[1];
rz(-2.5271509) q[1];
sx q[1];
rz(2.9628672) q[1];
rz(-pi) q[2];
rz(-0.68671642) q[3];
sx q[3];
rz(-2.2855008) q[3];
sx q[3];
rz(-1.2031632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19053261) q[2];
sx q[2];
rz(-0.6742) q[2];
sx q[2];
rz(2.8437957) q[2];
rz(2.8262302) q[3];
sx q[3];
rz(-0.30839977) q[3];
sx q[3];
rz(-1.4006008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067658871) q[0];
sx q[0];
rz(-0.44026259) q[0];
sx q[0];
rz(-0.0010781188) q[0];
rz(-3.0792117) q[1];
sx q[1];
rz(-1.1515113) q[1];
sx q[1];
rz(1.6463702) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1414215) q[0];
sx q[0];
rz(-2.6323491) q[0];
sx q[0];
rz(-1.4285136) q[0];
rz(-pi) q[1];
rz(-2.5220715) q[2];
sx q[2];
rz(-2.4378901) q[2];
sx q[2];
rz(1.0985398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8797258) q[1];
sx q[1];
rz(-0.84964439) q[1];
sx q[1];
rz(-2.7489064) q[1];
rz(-2.2132164) q[3];
sx q[3];
rz(-1.7084261) q[3];
sx q[3];
rz(-0.22882593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8794787) q[2];
sx q[2];
rz(-2.414119) q[2];
sx q[2];
rz(-1.0564085) q[2];
rz(-2.2004755) q[3];
sx q[3];
rz(-1.859715) q[3];
sx q[3];
rz(-0.077805422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62517405) q[0];
sx q[0];
rz(-0.70676261) q[0];
sx q[0];
rz(-1.9551552) q[0];
rz(0.9390074) q[1];
sx q[1];
rz(-0.72827315) q[1];
sx q[1];
rz(2.3659288) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.01919) q[0];
sx q[0];
rz(-1.3870326) q[0];
sx q[0];
rz(1.8792436) q[0];
rz(-pi) q[1];
rz(-0.93251419) q[2];
sx q[2];
rz(-2.4862977) q[2];
sx q[2];
rz(0.95791657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6104148) q[1];
sx q[1];
rz(-1.9912331) q[1];
sx q[1];
rz(-1.1320482) q[1];
rz(-pi) q[2];
rz(-1.6157563) q[3];
sx q[3];
rz(-1.3507843) q[3];
sx q[3];
rz(2.3156692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0607854) q[2];
sx q[2];
rz(-0.60254097) q[2];
sx q[2];
rz(1.5134643) q[2];
rz(1.8166517) q[3];
sx q[3];
rz(-2.5038268) q[3];
sx q[3];
rz(-0.12346867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.093928) q[0];
sx q[0];
rz(-1.0075189) q[0];
sx q[0];
rz(0.63543332) q[0];
rz(1.1134998) q[1];
sx q[1];
rz(-2.7058388) q[1];
sx q[1];
rz(-0.39438549) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6187103) q[0];
sx q[0];
rz(-3.0268207) q[0];
sx q[0];
rz(0.7392493) q[0];
rz(-pi) q[1];
rz(2.7761781) q[2];
sx q[2];
rz(-2.6435249) q[2];
sx q[2];
rz(0.77332449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0985165) q[1];
sx q[1];
rz(-1.9253007) q[1];
sx q[1];
rz(2.6621044) q[1];
x q[2];
rz(0.66310574) q[3];
sx q[3];
rz(-1.6702456) q[3];
sx q[3];
rz(-0.15583043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34027943) q[2];
sx q[2];
rz(-1.7873414) q[2];
sx q[2];
rz(-0.24564965) q[2];
rz(-0.73686016) q[3];
sx q[3];
rz(-2.1818678) q[3];
sx q[3];
rz(2.5514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067829475) q[0];
sx q[0];
rz(-1.2968061) q[0];
sx q[0];
rz(-0.8514362) q[0];
rz(-0.40359452) q[1];
sx q[1];
rz(-0.60786301) q[1];
sx q[1];
rz(3.0013989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64463413) q[0];
sx q[0];
rz(-1.2063364) q[0];
sx q[0];
rz(-0.24043997) q[0];
x q[1];
rz(-2.9459566) q[2];
sx q[2];
rz(-1.6802898) q[2];
sx q[2];
rz(-2.1957731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.043958943) q[1];
sx q[1];
rz(-1.2737899) q[1];
sx q[1];
rz(-0.87451248) q[1];
x q[2];
rz(2.5032205) q[3];
sx q[3];
rz(-2.4536096) q[3];
sx q[3];
rz(1.0561266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95789528) q[2];
sx q[2];
rz(-0.87527466) q[2];
sx q[2];
rz(0.20269205) q[2];
rz(-1.1787339) q[3];
sx q[3];
rz(-0.26114994) q[3];
sx q[3];
rz(1.096426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.390585) q[0];
sx q[0];
rz(-0.035476606) q[0];
sx q[0];
rz(-1.1998727) q[0];
rz(-1.9203145) q[1];
sx q[1];
rz(-0.55322909) q[1];
sx q[1];
rz(-1.663307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2182491) q[0];
sx q[0];
rz(-1.5594421) q[0];
sx q[0];
rz(1.5664805) q[0];
rz(0.37931077) q[2];
sx q[2];
rz(-2.0294445) q[2];
sx q[2];
rz(-0.5543602) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85046613) q[1];
sx q[1];
rz(-1.0134103) q[1];
sx q[1];
rz(-1.7859573) q[1];
rz(-pi) q[2];
rz(-2.8536232) q[3];
sx q[3];
rz(-1.4060296) q[3];
sx q[3];
rz(-2.3138177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.025909802) q[2];
sx q[2];
rz(-2.2175711) q[2];
sx q[2];
rz(2.8795418) q[2];
rz(2.986749) q[3];
sx q[3];
rz(-0.20379977) q[3];
sx q[3];
rz(-0.039732963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994951) q[0];
sx q[0];
rz(-2.4423548) q[0];
sx q[0];
rz(-2.5788838) q[0];
rz(-1.5744677) q[1];
sx q[1];
rz(-1.8561615) q[1];
sx q[1];
rz(-0.80975103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375423) q[0];
sx q[0];
rz(-1.3464377) q[0];
sx q[0];
rz(1.9363008) q[0];
rz(2.3793567) q[2];
sx q[2];
rz(-2.3272772) q[2];
sx q[2];
rz(-1.7310614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6470064) q[1];
sx q[1];
rz(-1.9368163) q[1];
sx q[1];
rz(-0.33509071) q[1];
x q[2];
rz(3.0122981) q[3];
sx q[3];
rz(-0.92513771) q[3];
sx q[3];
rz(-0.66619841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8030545) q[2];
sx q[2];
rz(-0.43236098) q[2];
sx q[2];
rz(-1.4114678) q[2];
rz(-0.159053) q[3];
sx q[3];
rz(-1.8820857) q[3];
sx q[3];
rz(2.4771396) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14801046) q[0];
sx q[0];
rz(-2.9382886) q[0];
sx q[0];
rz(-2.5116442) q[0];
rz(-2.3155164) q[1];
sx q[1];
rz(-0.45336205) q[1];
sx q[1];
rz(-1.0005383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2824018) q[0];
sx q[0];
rz(-2.1127208) q[0];
sx q[0];
rz(-1.6247092) q[0];
rz(-1.8597322) q[2];
sx q[2];
rz(-1.1232159) q[2];
sx q[2];
rz(-0.96842365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6521451) q[1];
sx q[1];
rz(-1.9838732) q[1];
sx q[1];
rz(2.4484171) q[1];
x q[2];
rz(-0.54686336) q[3];
sx q[3];
rz(-0.78843964) q[3];
sx q[3];
rz(2.2347286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91653812) q[2];
sx q[2];
rz(-0.30320898) q[2];
sx q[2];
rz(1.0603325) q[2];
rz(-0.66925085) q[3];
sx q[3];
rz(-2.4762912) q[3];
sx q[3];
rz(-0.76702142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6172248) q[0];
sx q[0];
rz(-1.1001294) q[0];
sx q[0];
rz(-0.97468162) q[0];
rz(-0.928448) q[1];
sx q[1];
rz(-1.4310373) q[1];
sx q[1];
rz(1.9824082) q[1];
rz(0.84865271) q[2];
sx q[2];
rz(-2.3046222) q[2];
sx q[2];
rz(0.23474856) q[2];
rz(-1.8228788) q[3];
sx q[3];
rz(-0.85267259) q[3];
sx q[3];
rz(0.55894077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
