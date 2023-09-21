OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(2.9705272) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(-2.7984518) q[1];
sx q[1];
rz(-1.8106102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(1.5050423) q[0];
rz(-pi) q[1];
rz(-1.5491484) q[2];
sx q[2];
rz(-1.0550371) q[2];
sx q[2];
rz(1.5390918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6273856) q[1];
sx q[1];
rz(-1.5537019) q[1];
sx q[1];
rz(-1.8219201) q[1];
x q[2];
rz(-1.347723) q[3];
sx q[3];
rz(-2.2017456) q[3];
sx q[3];
rz(2.0032721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.3295056) q[2];
rz(-1.7261516) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4703003) q[0];
sx q[0];
rz(-1.6912582) q[0];
sx q[0];
rz(2.5380773) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8882206) q[2];
sx q[2];
rz(-2.1183287) q[2];
sx q[2];
rz(0.65371338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2806432) q[1];
sx q[1];
rz(-1.3480535) q[1];
sx q[1];
rz(1.2692979) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2829451) q[3];
sx q[3];
rz(-0.24752366) q[3];
sx q[3];
rz(-1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(-2.7712908) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-2.8252576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4838594) q[0];
sx q[0];
rz(-1.3257926) q[0];
sx q[0];
rz(0.0018951456) q[0];
rz(2.5691368) q[2];
sx q[2];
rz(-1.9334031) q[2];
sx q[2];
rz(-0.87871694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0132732) q[1];
sx q[1];
rz(-1.1228021) q[1];
sx q[1];
rz(1.1303348) q[1];
x q[2];
rz(-2.8510677) q[3];
sx q[3];
rz(-1.2447262) q[3];
sx q[3];
rz(-2.3323004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8320223) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(-3.0679697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.081110031) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(0.62603355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44805376) q[1];
sx q[1];
rz(-1.5876019) q[1];
sx q[1];
rz(-2.0190713) q[1];
x q[2];
rz(-1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058218) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(3.0900893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52841016) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(0.86092237) q[0];
rz(-pi) q[1];
rz(-2.0578458) q[2];
sx q[2];
rz(-2.8205928) q[2];
sx q[2];
rz(3.0082862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67700779) q[1];
sx q[1];
rz(-1.1255463) q[1];
sx q[1];
rz(2.9693309) q[1];
x q[2];
rz(0.73370917) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(2.0165781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705868) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88887963) q[0];
sx q[0];
rz(-1.7656592) q[0];
sx q[0];
rz(1.3937508) q[0];
rz(-pi) q[1];
rz(-0.27481885) q[2];
sx q[2];
rz(-1.8289369) q[2];
sx q[2];
rz(-1.9354265) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6492918) q[1];
sx q[1];
rz(-1.2052844) q[1];
sx q[1];
rz(0.11458061) q[1];
rz(-pi) q[2];
rz(-2.5842651) q[3];
sx q[3];
rz(-1.3152272) q[3];
sx q[3];
rz(2.717358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.3827682) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7145342) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9761336) q[0];
sx q[0];
rz(-1.4384067) q[0];
sx q[0];
rz(-2.4112941) q[0];
rz(-pi) q[1];
rz(2.2824077) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(2.8189916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.84278216) q[1];
sx q[1];
rz(-2.543445) q[1];
sx q[1];
rz(0.36069718) q[1];
x q[2];
rz(-1.2956804) q[3];
sx q[3];
rz(-2.0886305) q[3];
sx q[3];
rz(0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-0.99532551) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156408) q[0];
sx q[0];
rz(-2.3693187) q[0];
sx q[0];
rz(0.036235972) q[0];
rz(-2.4023513) q[2];
sx q[2];
rz(-2.2695702) q[2];
sx q[2];
rz(-0.07585635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98098552) q[1];
sx q[1];
rz(-1.3541344) q[1];
sx q[1];
rz(0.76088455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4790672) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(-2.0446387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(2.8213815) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8378728) q[0];
sx q[0];
rz(-1.2178197) q[0];
sx q[0];
rz(-0.2210125) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9270883) q[2];
sx q[2];
rz(-2.5062343) q[2];
sx q[2];
rz(-0.84745849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.853211) q[1];
sx q[1];
rz(-0.92331159) q[1];
sx q[1];
rz(2.7700469) q[1];
rz(-pi) q[2];
rz(0.50001486) q[3];
sx q[3];
rz(-0.68272907) q[3];
sx q[3];
rz(2.6422215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(1.43191) q[0];
rz(2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(2.3419103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(0.090263788) q[0];
x q[1];
rz(1.6522371) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(2.75027) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4950072) q[3];
sx q[3];
rz(-1.4255376) q[3];
sx q[3];
rz(-0.69243542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(-2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-1.5699301) q[2];
sx q[2];
rz(-1.9067087) q[2];
sx q[2];
rz(-2.9839347) q[2];
rz(1.7351032) q[3];
sx q[3];
rz(-2.3621029) q[3];
sx q[3];
rz(2.7030871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];