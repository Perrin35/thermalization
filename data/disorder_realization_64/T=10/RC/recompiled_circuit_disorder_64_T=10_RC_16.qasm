OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927147) q[0];
sx q[0];
rz(-1.9034916) q[0];
sx q[0];
rz(-0.96941745) q[0];
rz(-1.9900436) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-0.56564769) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3803346) q[1];
sx q[1];
rz(-2.3654656) q[1];
sx q[1];
rz(-0.18591979) q[1];
x q[2];
rz(2.1943135) q[3];
sx q[3];
rz(-1.2926658) q[3];
sx q[3];
rz(2.5840685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-0.0016454776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5434108) q[0];
sx q[0];
rz(-3.0953005) q[0];
sx q[0];
rz(1.8020736) q[0];
x q[1];
rz(2.8029289) q[2];
sx q[2];
rz(-2.4709457) q[2];
sx q[2];
rz(1.8530958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.590608) q[1];
sx q[1];
rz(-1.814517) q[1];
sx q[1];
rz(2.3477712) q[1];
rz(1.315829) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(-0.69873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-2.3143342) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-2.0522096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339061) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-0.49467996) q[0];
rz(-pi) q[1];
x q[1];
rz(2.543534) q[2];
sx q[2];
rz(-1.2130249) q[2];
sx q[2];
rz(-1.2146815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8100064) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(-1.3473131) q[1];
x q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.9867992) q[3];
sx q[3];
rz(-2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(0.35983905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4041876) q[0];
sx q[0];
rz(-2.2149937) q[0];
sx q[0];
rz(-0.23097158) q[0];
x q[1];
rz(1.4570518) q[2];
sx q[2];
rz(-2.2188088) q[2];
sx q[2];
rz(1.5485473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6872014) q[1];
sx q[1];
rz(-1.9481716) q[1];
sx q[1];
rz(-2.6881933) q[1];
rz(-0.30712819) q[3];
sx q[3];
rz(-1.9046475) q[3];
sx q[3];
rz(3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9969479) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(1.0353154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0400378) q[0];
sx q[0];
rz(-0.77417513) q[0];
sx q[0];
rz(0.9206307) q[0];
rz(1.5315227) q[2];
sx q[2];
rz(-1.0796667) q[2];
sx q[2];
rz(-0.66205762) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94381911) q[1];
sx q[1];
rz(-0.10995956) q[1];
sx q[1];
rz(-2.6206559) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4953793) q[3];
sx q[3];
rz(-2.8013902) q[3];
sx q[3];
rz(-0.37341213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(2.8055577) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(-1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(0.2072269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6672872) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(-2.5179203) q[0];
x q[1];
rz(0.60116641) q[2];
sx q[2];
rz(-1.9631557) q[2];
sx q[2];
rz(-0.62077921) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.065437) q[1];
sx q[1];
rz(-1.3216615) q[1];
sx q[1];
rz(2.8912828) q[1];
rz(-pi) q[2];
rz(-0.16222555) q[3];
sx q[3];
rz(-2.3174006) q[3];
sx q[3];
rz(1.1713067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(0.2917372) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9852108) q[0];
sx q[0];
rz(-1.1522326) q[0];
sx q[0];
rz(1.5147989) q[0];
rz(-pi) q[1];
rz(-1.8740011) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(2.6333957) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8818672) q[1];
sx q[1];
rz(-1.423466) q[1];
sx q[1];
rz(2.8357361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6187906) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-2.4777381) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(-2.1829139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6129235) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(0.59707609) q[0];
rz(1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(1.4657071) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9194591) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(1.7807351) q[1];
x q[2];
rz(2.9577191) q[3];
sx q[3];
rz(-0.5920147) q[3];
sx q[3];
rz(-1.5541935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0083996) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(-1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(0.46554309) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(2.4337016) q[0];
x q[1];
rz(1.8579673) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(0.62921333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8852639) q[1];
sx q[1];
rz(-1.1210821) q[1];
sx q[1];
rz(1.9535669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2113308) q[3];
sx q[3];
rz(-1.8065479) q[3];
sx q[3];
rz(-0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(2.8695316) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96888992) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(1.885528) q[0];
x q[1];
rz(1.2920612) q[2];
sx q[2];
rz(-0.50694743) q[2];
sx q[2];
rz(2.7915733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-2.1464159) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9243745) q[3];
sx q[3];
rz(-1.5434885) q[3];
sx q[3];
rz(-2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5738752) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(0.81417685) q[2];
sx q[2];
rz(-1.3635175) q[2];
sx q[2];
rz(-1.2843532) q[2];
rz(2.3253757) q[3];
sx q[3];
rz(-2.4945138) q[3];
sx q[3];
rz(2.1993779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
