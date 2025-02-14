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
rz(-2.8060198) q[0];
sx q[0];
rz(-1.6874474) q[0];
sx q[0];
rz(2.1079221) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.570734) q[0];
sx q[0];
rz(-1.6964127) q[0];
sx q[0];
rz(-1.7616338) q[0];
rz(-pi) q[1];
rz(-1.6160585) q[2];
sx q[2];
rz(-3.0590364) q[2];
sx q[2];
rz(2.1545067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42130175) q[1];
sx q[1];
rz(-0.20245782) q[1];
sx q[1];
rz(2.5295707) q[1];
rz(3.1107424) q[3];
sx q[3];
rz(-2.9992963) q[3];
sx q[3];
rz(0.22103413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7734163) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(0.82484335) q[2];
rz(-2.8090737) q[3];
sx q[3];
rz(-2.2907292) q[3];
sx q[3];
rz(-0.48377812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9124209) q[0];
sx q[0];
rz(-3.0569172) q[0];
sx q[0];
rz(-0.53572768) q[0];
rz(-2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.7492693) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014382) q[0];
sx q[0];
rz(-1.2291698) q[0];
sx q[0];
rz(-1.0761976) q[0];
rz(0.34113348) q[2];
sx q[2];
rz(-1.5798414) q[2];
sx q[2];
rz(1.1760528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0619095) q[1];
sx q[1];
rz(-2.3259729) q[1];
sx q[1];
rz(-3.1119926) q[1];
rz(-2.4588938) q[3];
sx q[3];
rz(-0.16541741) q[3];
sx q[3];
rz(1.1755652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2367737) q[2];
sx q[2];
rz(-2.3084013) q[2];
sx q[2];
rz(2.0429677) q[2];
rz(-2.3084579) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5611834) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(-0.68823254) q[0];
rz(1.2511823) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(-1.2122663) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466169) q[0];
sx q[0];
rz(-2.2612834) q[0];
sx q[0];
rz(0.63027905) q[0];
x q[1];
rz(-1.8463676) q[2];
sx q[2];
rz(-1.3418806) q[2];
sx q[2];
rz(0.51425291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7292068) q[1];
sx q[1];
rz(-2.4826035) q[1];
sx q[1];
rz(3.1324196) q[1];
rz(-1.6501649) q[3];
sx q[3];
rz(-1.575633) q[3];
sx q[3];
rz(-1.2629379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7324098) q[2];
sx q[2];
rz(-2.1467291) q[2];
sx q[2];
rz(0.44962064) q[2];
rz(-0.85363394) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(-0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39321008) q[0];
sx q[0];
rz(-1.0164096) q[0];
sx q[0];
rz(-2.7384695) q[0];
rz(-2.368811) q[1];
sx q[1];
rz(-1.2330331) q[1];
sx q[1];
rz(-0.82537878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8050418) q[0];
sx q[0];
rz(-0.8835578) q[0];
sx q[0];
rz(-0.42091333) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4140777) q[2];
sx q[2];
rz(-0.036002654) q[2];
sx q[2];
rz(1.9787479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.014848901) q[1];
sx q[1];
rz(-2.150661) q[1];
sx q[1];
rz(0.090437263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.230497) q[3];
sx q[3];
rz(-2.0269217) q[3];
sx q[3];
rz(-0.62247372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1021425) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(0.5086745) q[2];
rz(1.8968808) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(-1.8083474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4077069) q[0];
sx q[0];
rz(-1.5579959) q[0];
sx q[0];
rz(2.7852614) q[0];
rz(-1.412926) q[1];
sx q[1];
rz(-1.8316385) q[1];
sx q[1];
rz(1.8467356) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0771478) q[0];
sx q[0];
rz(-1.1379717) q[0];
sx q[0];
rz(-1.8124142) q[0];
rz(-0.6614954) q[2];
sx q[2];
rz(-1.7057888) q[2];
sx q[2];
rz(2.1141426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2834245) q[1];
sx q[1];
rz(-1.3407665) q[1];
sx q[1];
rz(2.7344804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44900198) q[3];
sx q[3];
rz(-2.0709566) q[3];
sx q[3];
rz(-1.4003818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20092189) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(2.5015639) q[2];
rz(-0.43404964) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(-2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(2.3906999) q[0];
rz(2.3789876) q[1];
sx q[1];
rz(-1.4354939) q[1];
sx q[1];
rz(-2.6079752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92350875) q[0];
sx q[0];
rz(-1.8663915) q[0];
sx q[0];
rz(-1.8718998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0178823) q[2];
sx q[2];
rz(-2.3229282) q[2];
sx q[2];
rz(-2.1731203) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4590964) q[1];
sx q[1];
rz(-2.4630416) q[1];
sx q[1];
rz(1.7595923) q[1];
x q[2];
rz(0.98762767) q[3];
sx q[3];
rz(-2.0500053) q[3];
sx q[3];
rz(-0.85239172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0642455) q[2];
sx q[2];
rz(-0.40234819) q[2];
sx q[2];
rz(2.2557491) q[2];
rz(1.3028076) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9293905) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-2.5481664) q[0];
rz(-1.628283) q[1];
sx q[1];
rz(-1.1791041) q[1];
sx q[1];
rz(-2.5679307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79686576) q[0];
sx q[0];
rz(-1.0930976) q[0];
sx q[0];
rz(1.9928689) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3232347) q[2];
sx q[2];
rz(-2.1344961) q[2];
sx q[2];
rz(-2.7147788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57854125) q[1];
sx q[1];
rz(-2.6669901) q[1];
sx q[1];
rz(2.5792349) q[1];
rz(1.6214293) q[3];
sx q[3];
rz(-1.4922499) q[3];
sx q[3];
rz(0.95548624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71184984) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(1.708606) q[2];
rz(-1.9728707) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(1.6247113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.8716005) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(0.20371833) q[0];
rz(-1.8261955) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(-1.3320097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34623554) q[0];
sx q[0];
rz(-0.95690823) q[0];
sx q[0];
rz(2.3625663) q[0];
rz(-pi) q[1];
rz(0.90317995) q[2];
sx q[2];
rz(-2.0364663) q[2];
sx q[2];
rz(-2.6671034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9086919) q[1];
sx q[1];
rz(-2.3756643) q[1];
sx q[1];
rz(-0.58118622) q[1];
rz(-0.12961403) q[3];
sx q[3];
rz(-1.4091421) q[3];
sx q[3];
rz(1.1974789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
sx q[2];
rz(1.7601684) q[2];
rz(2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(-0.74530017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7119174) q[0];
sx q[0];
rz(-1.6918007) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(-1.3445541) q[1];
sx q[1];
rz(-2.5090736) q[1];
sx q[1];
rz(1.9154027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525402) q[0];
sx q[0];
rz(-1.2864094) q[0];
sx q[0];
rz(-3.000598) q[0];
rz(2.7550952) q[2];
sx q[2];
rz(-1.1697239) q[2];
sx q[2];
rz(-2.3569466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74599712) q[1];
sx q[1];
rz(-1.5052553) q[1];
sx q[1];
rz(2.3941674) q[1];
x q[2];
rz(0.33554797) q[3];
sx q[3];
rz(-1.0764383) q[3];
sx q[3];
rz(1.7691607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58174497) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(-1.6017412) q[2];
rz(1.7848484) q[3];
sx q[3];
rz(-2.609085) q[3];
sx q[3];
rz(1.9073568) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8643879) q[0];
sx q[0];
rz(-1.5978403) q[0];
sx q[0];
rz(-0.10449115) q[0];
rz(-1.3970207) q[1];
sx q[1];
rz(-1.7897768) q[1];
sx q[1];
rz(-1.6207961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6817271) q[0];
sx q[0];
rz(-2.2978373) q[0];
sx q[0];
rz(-1.8521502) q[0];
x q[1];
rz(-2.7149523) q[2];
sx q[2];
rz(-0.37084118) q[2];
sx q[2];
rz(-1.1877354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11776671) q[1];
sx q[1];
rz(-1.858288) q[1];
sx q[1];
rz(-1.6790798) q[1];
rz(3.0836068) q[3];
sx q[3];
rz(-2.1515905) q[3];
sx q[3];
rz(0.58422663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4349159) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(1.288877) q[2];
rz(-1.0051109) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(-3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91916753) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(-2.2617321) q[1];
sx q[1];
rz(-2.6251371) q[1];
sx q[1];
rz(0.62216204) q[1];
rz(0.11779412) q[2];
sx q[2];
rz(-0.5024903) q[2];
sx q[2];
rz(2.6033664) q[2];
rz(-2.5414657) q[3];
sx q[3];
rz(-1.5059581) q[3];
sx q[3];
rz(-3.0653421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
