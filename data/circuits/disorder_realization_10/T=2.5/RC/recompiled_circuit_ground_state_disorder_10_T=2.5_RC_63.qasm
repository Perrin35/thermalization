OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(-0.9932819) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(-1.7239404) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7844789) q[0];
sx q[0];
rz(-2.3001101) q[0];
sx q[0];
rz(-0.48805663) q[0];
rz(-0.63739802) q[2];
sx q[2];
rz(-2.3294724) q[2];
sx q[2];
rz(-0.82773436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65517814) q[1];
sx q[1];
rz(-1.8433807) q[1];
sx q[1];
rz(0.64719871) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7050691) q[3];
sx q[3];
rz(-1.6032919) q[3];
sx q[3];
rz(2.1509347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.950497) q[2];
sx q[2];
rz(-1.82093) q[2];
sx q[2];
rz(2.0868059) q[2];
rz(0.47109207) q[3];
sx q[3];
rz(-1.4548917) q[3];
sx q[3];
rz(-2.2733222) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3926113) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(2.9066322) q[0];
rz(1.2745534) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(0.68769208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4859568) q[0];
sx q[0];
rz(-1.63795) q[0];
sx q[0];
rz(-0.32819076) q[0];
rz(0.85854395) q[2];
sx q[2];
rz(-1.9496103) q[2];
sx q[2];
rz(1.0114947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60581698) q[1];
sx q[1];
rz(-1.5652085) q[1];
sx q[1];
rz(1.5833012) q[1];
x q[2];
rz(-3.0491203) q[3];
sx q[3];
rz(-2.1824129) q[3];
sx q[3];
rz(-0.05449748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4253652) q[2];
sx q[2];
rz(-1.3024412) q[2];
sx q[2];
rz(-2.9322374) q[2];
rz(-0.9730722) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.3000325) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(5/(11*pi)) q[0];
rz(1.2380098) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(-0.091781052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9301355) q[0];
sx q[0];
rz(-1.2263023) q[0];
sx q[0];
rz(-2.1449367) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0726542) q[2];
sx q[2];
rz(-1.6668011) q[2];
sx q[2];
rz(-1.7092741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8615177) q[1];
sx q[1];
rz(-2.8535378) q[1];
sx q[1];
rz(2.3029598) q[1];
rz(-pi) q[2];
rz(1.1385659) q[3];
sx q[3];
rz(-1.5097268) q[3];
sx q[3];
rz(-1.40738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38516513) q[2];
sx q[2];
rz(-1.5084718) q[2];
sx q[2];
rz(-2.7623994) q[2];
rz(2.3181629) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(-0.2894952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32651153) q[0];
sx q[0];
rz(-0.87711763) q[0];
sx q[0];
rz(-0.99863482) q[0];
rz(-0.48209349) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(1.3179717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362396) q[0];
sx q[0];
rz(-2.0990206) q[0];
sx q[0];
rz(1.1861237) q[0];
rz(-0.58546328) q[2];
sx q[2];
rz(-1.778086) q[2];
sx q[2];
rz(-0.49400615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6819344) q[1];
sx q[1];
rz(-1.0409969) q[1];
sx q[1];
rz(1.3270686) q[1];
rz(-2.34818) q[3];
sx q[3];
rz(-1.9559091) q[3];
sx q[3];
rz(-2.957893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92129293) q[2];
sx q[2];
rz(-0.90196323) q[2];
sx q[2];
rz(-0.48545066) q[2];
rz(-0.38393936) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(0.78401047) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081472814) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(-0.01165788) q[0];
rz(3.0043789) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(1.8395909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6377651) q[0];
sx q[0];
rz(-1.3877467) q[0];
sx q[0];
rz(1.6812912) q[0];
rz(-pi) q[1];
rz(-2.2574591) q[2];
sx q[2];
rz(-2.4120286) q[2];
sx q[2];
rz(-1.8232249) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3350388) q[1];
sx q[1];
rz(-1.8825899) q[1];
sx q[1];
rz(0.11374082) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7632159) q[3];
sx q[3];
rz(-1.9073745) q[3];
sx q[3];
rz(-0.61016309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4181218) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(-2.8483025) q[2];
rz(-3.0784741) q[3];
sx q[3];
rz(-1.9189574) q[3];
sx q[3];
rz(2.2466834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191331) q[0];
sx q[0];
rz(-2.8526511) q[0];
sx q[0];
rz(3.1030848) q[0];
rz(-1.0844082) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(0.96673059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.777323) q[0];
sx q[0];
rz(-1.7614953) q[0];
sx q[0];
rz(2.9421115) q[0];
x q[1];
rz(1.6183774) q[2];
sx q[2];
rz(-0.26740012) q[2];
sx q[2];
rz(-2.9427955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37361426) q[1];
sx q[1];
rz(-1.9100338) q[1];
sx q[1];
rz(-1.148073) q[1];
x q[2];
rz(2.8916675) q[3];
sx q[3];
rz(-0.93533134) q[3];
sx q[3];
rz(-0.15033406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4806369) q[2];
sx q[2];
rz(-2.9065242) q[2];
sx q[2];
rz(0.98141518) q[2];
rz(-1.6437982) q[3];
sx q[3];
rz(-1.2621597) q[3];
sx q[3];
rz(1.5463382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.3129231) q[0];
sx q[0];
rz(-2.4607615) q[0];
sx q[0];
rz(2.8269826) q[0];
rz(-2.9529052) q[1];
sx q[1];
rz(-2.7068832) q[1];
sx q[1];
rz(-0.46868086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.311694) q[0];
sx q[0];
rz(-1.7009987) q[0];
sx q[0];
rz(-0.24451406) q[0];
rz(-3.1397083) q[2];
sx q[2];
rz(-1.7385637) q[2];
sx q[2];
rz(1.1275856) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3568253) q[1];
sx q[1];
rz(-1.2677437) q[1];
sx q[1];
rz(0.81929889) q[1];
x q[2];
rz(-1.8953034) q[3];
sx q[3];
rz(-1.1238928) q[3];
sx q[3];
rz(2.383142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9357052) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(0.97071281) q[2];
rz(-2.6054221) q[3];
sx q[3];
rz(-1.9325117) q[3];
sx q[3];
rz(0.71603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69990528) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(-1.6957138) q[0];
rz(2.1031759) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(-0.081092484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2807407) q[0];
sx q[0];
rz(-2.3978455) q[0];
sx q[0];
rz(-2.70983) q[0];
x q[1];
rz(-1.2732701) q[2];
sx q[2];
rz(-0.93704455) q[2];
sx q[2];
rz(0.81281205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7837569) q[1];
sx q[1];
rz(-1.862613) q[1];
sx q[1];
rz(3.0213657) q[1];
x q[2];
rz(2.0399953) q[3];
sx q[3];
rz(-2.0939671) q[3];
sx q[3];
rz(1.8947471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1991835) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(-0.24277631) q[2];
rz(-1.5593922) q[3];
sx q[3];
rz(-2.715761) q[3];
sx q[3];
rz(-1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9026069) q[0];
sx q[0];
rz(-2.1098397) q[0];
sx q[0];
rz(1.8865939) q[0];
rz(-1.7651419) q[1];
sx q[1];
rz(-1.0245198) q[1];
sx q[1];
rz(0.66744101) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7783035) q[0];
sx q[0];
rz(-0.61909715) q[0];
sx q[0];
rz(0.97514345) q[0];
x q[1];
rz(-3.0632067) q[2];
sx q[2];
rz(-2.6437685) q[2];
sx q[2];
rz(-2.5428685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.046958663) q[1];
sx q[1];
rz(-0.65494767) q[1];
sx q[1];
rz(1.0818693) q[1];
rz(-pi) q[2];
rz(0.38102229) q[3];
sx q[3];
rz(-1.8748611) q[3];
sx q[3];
rz(-1.1883433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5218375) q[2];
sx q[2];
rz(-0.73885584) q[2];
sx q[2];
rz(0.45836207) q[2];
rz(-1.5357664) q[3];
sx q[3];
rz(-1.0777487) q[3];
sx q[3];
rz(2.2569236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.89852029) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(0.36548734) q[0];
rz(-0.99994031) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(1.4568636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3709049) q[0];
sx q[0];
rz(-1.7720776) q[0];
sx q[0];
rz(-2.1343436) q[0];
x q[1];
rz(-1.959895) q[2];
sx q[2];
rz(-0.28841296) q[2];
sx q[2];
rz(2.4450977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36023203) q[1];
sx q[1];
rz(-1.39729) q[1];
sx q[1];
rz(0.75299112) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6014464) q[3];
sx q[3];
rz(-2.3688101) q[3];
sx q[3];
rz(-2.5908822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(-2.136039) q[2];
rz(-2.4908861) q[3];
sx q[3];
rz(-1.4395827) q[3];
sx q[3];
rz(-2.2910291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093813048) q[0];
sx q[0];
rz(-2.35738) q[0];
sx q[0];
rz(2.1398075) q[0];
rz(-1.7306937) q[1];
sx q[1];
rz(-1.4374562) q[1];
sx q[1];
rz(1.1387574) q[1];
rz(-1.2415573) q[2];
sx q[2];
rz(-1.4598979) q[2];
sx q[2];
rz(2.0171063) q[2];
rz(-3.0372314) q[3];
sx q[3];
rz(-1.6634533) q[3];
sx q[3];
rz(-2.8692393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
