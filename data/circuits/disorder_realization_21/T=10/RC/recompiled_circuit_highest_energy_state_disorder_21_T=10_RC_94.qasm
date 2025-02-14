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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(2.5330438) q[0];
rz(-0.67148709) q[1];
sx q[1];
rz(3.804764) q[1];
sx q[1];
rz(10.884203) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08069399) q[0];
sx q[0];
rz(-2.0804724) q[0];
sx q[0];
rz(-0.26943785) q[0];
rz(-2.0107277) q[2];
sx q[2];
rz(-0.49060433) q[2];
sx q[2];
rz(-3.0705796) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8512177) q[1];
sx q[1];
rz(-0.77271116) q[1];
sx q[1];
rz(0.98224838) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75992775) q[3];
sx q[3];
rz(-1.6960521) q[3];
sx q[3];
rz(1.4533991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3964316) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(0.82610899) q[2];
rz(2.7769026) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(2.8472692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18225886) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(0.44573319) q[0];
rz(-0.590473) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(2.7139434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70156258) q[0];
sx q[0];
rz(-2.0272581) q[0];
sx q[0];
rz(0.052340074) q[0];
x q[1];
rz(2.580709) q[2];
sx q[2];
rz(-2.2429401) q[2];
sx q[2];
rz(2.667556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1688437) q[1];
sx q[1];
rz(-0.66509897) q[1];
sx q[1];
rz(1.14023) q[1];
rz(1.5656972) q[3];
sx q[3];
rz(-2.3672464) q[3];
sx q[3];
rz(2.4536595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.4599956) q[2];
sx q[2];
rz(-2.5619016) q[2];
sx q[2];
rz(-2.106529) q[2];
rz(-1.5486108) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(-0.88919324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11397938) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(-3.0300544) q[0];
rz(-1.5533718) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-2.9389971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3521393) q[0];
sx q[0];
rz(-1.7048262) q[0];
sx q[0];
rz(-1.171597) q[0];
rz(0.79988806) q[2];
sx q[2];
rz(-2.7383907) q[2];
sx q[2];
rz(0.099911913) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2059171) q[1];
sx q[1];
rz(-0.96921235) q[1];
sx q[1];
rz(2.4755073) q[1];
rz(-pi) q[2];
rz(2.5532365) q[3];
sx q[3];
rz(-0.6783456) q[3];
sx q[3];
rz(2.7916186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4488039) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(-2.1602574) q[2];
rz(2.8547309) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(-2.9441492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2364748) q[0];
sx q[0];
rz(-3.0906257) q[0];
sx q[0];
rz(0.68848759) q[0];
rz(2.9972637) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(2.3932638) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37667021) q[0];
sx q[0];
rz(-2.2895396) q[0];
sx q[0];
rz(-1.5115949) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2690721) q[2];
sx q[2];
rz(-1.3137443) q[2];
sx q[2];
rz(-1.014705) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1074368) q[1];
sx q[1];
rz(-1.7567051) q[1];
sx q[1];
rz(2.8380945) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9311936) q[3];
sx q[3];
rz(-2.4080417) q[3];
sx q[3];
rz(1.4195132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.051108483) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(-3.0221853) q[3];
sx q[3];
rz(-1.2179008) q[3];
sx q[3];
rz(2.8363805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3731821) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(0.61491948) q[0];
rz(2.3070295) q[1];
sx q[1];
rz(-1.395023) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1008719) q[0];
sx q[0];
rz(-0.03558579) q[0];
sx q[0];
rz(1.8338704) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70186285) q[2];
sx q[2];
rz(-0.80282513) q[2];
sx q[2];
rz(2.7868164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78042049) q[1];
sx q[1];
rz(-2.5675312) q[1];
sx q[1];
rz(1.3657544) q[1];
rz(0.74252601) q[3];
sx q[3];
rz(-2.1603185) q[3];
sx q[3];
rz(1.7537011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(0.92030805) q[2];
rz(-0.80095428) q[3];
sx q[3];
rz(-0.52217537) q[3];
sx q[3];
rz(-0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1717218) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(0.26611662) q[0];
rz(-1.4316106) q[1];
sx q[1];
rz(-1.1191198) q[1];
sx q[1];
rz(0.59763479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1475567) q[0];
sx q[0];
rz(-1.3809675) q[0];
sx q[0];
rz(0.61845431) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0805239) q[2];
sx q[2];
rz(-2.960223) q[2];
sx q[2];
rz(0.86038113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2955851) q[1];
sx q[1];
rz(-0.69900988) q[1];
sx q[1];
rz(3.0547804) q[1];
rz(-pi) q[2];
x q[2];
rz(1.131874) q[3];
sx q[3];
rz(-0.35177975) q[3];
sx q[3];
rz(-1.9466512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0133682) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(-2.5426148) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-2.7003761) q[3];
sx q[3];
rz(2.0986572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041158572) q[0];
sx q[0];
rz(-0.91223311) q[0];
sx q[0];
rz(-0.31901932) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.7496505) q[1];
sx q[1];
rz(-2.9999733) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7975811) q[0];
sx q[0];
rz(-2.251029) q[0];
sx q[0];
rz(1.4466831) q[0];
x q[1];
rz(0.0044195375) q[2];
sx q[2];
rz(-1.0194155) q[2];
sx q[2];
rz(-2.1194292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31894058) q[1];
sx q[1];
rz(-1.7590862) q[1];
sx q[1];
rz(0.55231185) q[1];
rz(-pi) q[2];
rz(1.0504248) q[3];
sx q[3];
rz(-1.4401575) q[3];
sx q[3];
rz(1.4121233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(-0.11036135) q[2];
rz(2.6335671) q[3];
sx q[3];
rz(-3.09943) q[3];
sx q[3];
rz(1.9917816) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-2.4899794) q[0];
sx q[0];
rz(-1.9223258) q[0];
rz(-3.1189175) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(-2.5882744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1071725) q[0];
sx q[0];
rz(-1.4822472) q[0];
sx q[0];
rz(2.2293985) q[0];
rz(-2.7560541) q[2];
sx q[2];
rz(-2.1111167) q[2];
sx q[2];
rz(-1.413942) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4412429) q[1];
sx q[1];
rz(-1.9158331) q[1];
sx q[1];
rz(2.5923968) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24711547) q[3];
sx q[3];
rz(-1.860642) q[3];
sx q[3];
rz(-0.71699504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(1.7919398) q[2];
rz(-3.0810629) q[3];
sx q[3];
rz(-2.2318201) q[3];
sx q[3];
rz(0.6282261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.26199207) q[0];
sx q[0];
rz(-0.99402004) q[0];
sx q[0];
rz(-0.0066198786) q[0];
rz(0.63893843) q[1];
sx q[1];
rz(-2.9220351) q[1];
sx q[1];
rz(-2.7117597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98243356) q[0];
sx q[0];
rz(-0.52839708) q[0];
sx q[0];
rz(2.3702752) q[0];
rz(1.1388135) q[2];
sx q[2];
rz(-1.8021934) q[2];
sx q[2];
rz(2.9803986) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62334261) q[1];
sx q[1];
rz(-2.3584705) q[1];
sx q[1];
rz(0.030695774) q[1];
rz(-pi) q[2];
rz(-2.5561629) q[3];
sx q[3];
rz(-1.6282035) q[3];
sx q[3];
rz(-1.2186183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.1755875) q[2];
sx q[2];
rz(3.1268934) q[2];
rz(-2.0059351) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(-1.3057825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.084918) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(1.2925451) q[0];
rz(2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(0.83052105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6097858) q[0];
sx q[0];
rz(-1.6702819) q[0];
sx q[0];
rz(-2.2983944) q[0];
x q[1];
rz(-1.8761329) q[2];
sx q[2];
rz(-2.3407901) q[2];
sx q[2];
rz(0.82894553) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78962425) q[1];
sx q[1];
rz(-1.3543173) q[1];
sx q[1];
rz(-3.0596872) q[1];
rz(-pi) q[2];
rz(-0.51955207) q[3];
sx q[3];
rz(-2.7249935) q[3];
sx q[3];
rz(2.8263457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(-2.9923934) q[3];
sx q[3];
rz(-2.7982893) q[3];
sx q[3];
rz(0.33299115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0014342) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(-0.99826605) q[1];
sx q[1];
rz(-1.2405735) q[1];
sx q[1];
rz(2.1539198) q[1];
rz(0.34299359) q[2];
sx q[2];
rz(-1.3806812) q[2];
sx q[2];
rz(1.5594202) q[2];
rz(-1.042749) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
