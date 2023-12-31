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
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9026731) q[0];
sx q[0];
rz(-2.9905149) q[0];
sx q[0];
rz(2.8047049) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2968729) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(-0.5288045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2831813) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(-1.7130501) q[1];
rz(-pi) q[2];
rz(1.7232056) q[3];
sx q[3];
rz(-1.8525575) q[3];
sx q[3];
rz(-1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.2715682) q[0];
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
rz(2.9257665) q[0];
x q[1];
rz(1.7550811) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(-2.266303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8010506) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1249173) q[3];
sx q[3];
rz(-1.8801873) q[3];
sx q[3];
rz(-0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19757195) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(2.5463085) q[0];
rz(1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-0.30644882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0551128) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(-2.8716645) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62927411) q[3];
sx q[3];
rz(-2.4089703) q[3];
sx q[3];
rz(2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(-2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(3.0217357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0474931) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(-3.0981307) q[0];
rz(2.7808431) q[2];
sx q[2];
rz(-1.7067688) q[2];
sx q[2];
rz(0.13775682) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3164697) q[1];
sx q[1];
rz(-2.9087062) q[1];
sx q[1];
rz(1.4941925) q[1];
x q[2];
rz(-0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64802058) q[0];
sx q[0];
rz(-2.9020712) q[0];
sx q[0];
rz(3.0859205) q[0];
rz(-2.901022) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(-2.0672928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8023194) q[1];
sx q[1];
rz(-2.5335651) q[1];
sx q[1];
rz(0.1165216) q[1];
rz(-pi) q[2];
rz(-0.48360444) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(-0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.038625) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(-1.6606768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(2.6631151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(-2.900219) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.092993) q[3];
sx q[3];
rz(-0.89266333) q[3];
sx q[3];
rz(-0.62852678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-2.1140153) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0425456) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(3.0142473) q[0];
rz(-pi) q[1];
rz(-0.84682805) q[2];
sx q[2];
rz(-0.70038751) q[2];
sx q[2];
rz(-2.2221153) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.160497) q[1];
x q[2];
rz(-2.2277101) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(-2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-2.1477551) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3562718) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(-2.450599) q[0];
x q[1];
rz(-0.71596594) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(1.3032608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.928927) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(-3.0287663) q[1];
rz(-0.21059489) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(-1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-2.0588493) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219014) q[0];
sx q[0];
rz(-1.1150517) q[0];
sx q[0];
rz(2.4569608) q[0];
rz(-pi) q[1];
rz(1.3547782) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(-0.87460364) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5987451) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(-0.41136841) q[1];
rz(-pi) q[2];
rz(2.2563124) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-1.0202927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-0.68181109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17446974) q[0];
sx q[0];
rz(-2.5645064) q[0];
sx q[0];
rz(1.8814927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85068662) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0360003) q[1];
sx q[1];
rz(-2.3899374) q[1];
sx q[1];
rz(-2.4113301) q[1];
rz(-pi) q[2];
rz(-0.52313157) q[3];
sx q[3];
rz(-1.0732068) q[3];
sx q[3];
rz(2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-2.5088359) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(-1.7726462) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
