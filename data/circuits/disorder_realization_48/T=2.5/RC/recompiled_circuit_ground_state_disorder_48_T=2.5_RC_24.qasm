OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(2.7383374) q[0];
rz(-0.23303214) q[1];
sx q[1];
rz(-1.4401399) q[1];
sx q[1];
rz(0.22388248) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163607) q[0];
sx q[0];
rz(-2.4852042) q[0];
sx q[0];
rz(-0.34968495) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2751132) q[2];
sx q[2];
rz(-1.472855) q[2];
sx q[2];
rz(2.3549454) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2157872) q[1];
sx q[1];
rz(-2.4443382) q[1];
sx q[1];
rz(1.9250735) q[1];
rz(-pi) q[2];
rz(2.4591962) q[3];
sx q[3];
rz(-0.93554893) q[3];
sx q[3];
rz(-1.531383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69748059) q[2];
sx q[2];
rz(-2.8287502) q[2];
sx q[2];
rz(0.035813896) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(-0.20619503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25074729) q[0];
sx q[0];
rz(-0.99531168) q[0];
sx q[0];
rz(-2.9440951) q[0];
rz(-0.47710553) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(-0.8055996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2700078) q[0];
sx q[0];
rz(-0.91579899) q[0];
sx q[0];
rz(-1.3359469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16519157) q[2];
sx q[2];
rz(-0.38536638) q[2];
sx q[2];
rz(-1.8866273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9979084) q[1];
sx q[1];
rz(-2.152118) q[1];
sx q[1];
rz(-0.79400009) q[1];
x q[2];
rz(-0.15112215) q[3];
sx q[3];
rz(-1.7100705) q[3];
sx q[3];
rz(2.3230011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3104559) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(-1.7384701) q[2];
rz(2.5095615) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(-0.12701756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8359351) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(-2.054731) q[0];
rz(2.944259) q[1];
sx q[1];
rz(-1.6355762) q[1];
sx q[1];
rz(0.33263439) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64651782) q[0];
sx q[0];
rz(-1.7103638) q[0];
sx q[0];
rz(1.1366751) q[0];
x q[1];
rz(-1.5015177) q[2];
sx q[2];
rz(-0.73613077) q[2];
sx q[2];
rz(2.2683805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.931529) q[1];
sx q[1];
rz(-1.3547661) q[1];
sx q[1];
rz(0.17976168) q[1];
rz(-pi) q[2];
rz(-0.72686355) q[3];
sx q[3];
rz(-1.3976025) q[3];
sx q[3];
rz(1.7693335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7670333) q[2];
sx q[2];
rz(-1.5908396) q[2];
sx q[2];
rz(1.4683051) q[2];
rz(0.0091920216) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(-2.3495242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62995768) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(-1.4427503) q[0];
rz(-2.1171782) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(0.82694298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31931811) q[0];
sx q[0];
rz(-1.3869023) q[0];
sx q[0];
rz(3.0660423) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2587772) q[2];
sx q[2];
rz(-1.4512296) q[2];
sx q[2];
rz(1.9012698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.045438529) q[1];
sx q[1];
rz(-1.3609582) q[1];
sx q[1];
rz(0.44053284) q[1];
rz(-2.2531281) q[3];
sx q[3];
rz(-1.4766221) q[3];
sx q[3];
rz(1.7674131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3247165) q[2];
sx q[2];
rz(-0.78410316) q[2];
sx q[2];
rz(-1.8640222) q[2];
rz(-0.68814021) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-0.96243206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14030305) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-2.9241614) q[0];
rz(-2.8561719) q[1];
sx q[1];
rz(-0.60573429) q[1];
sx q[1];
rz(-2.6434456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.065154) q[0];
sx q[0];
rz(-0.74238837) q[0];
sx q[0];
rz(-2.0641293) q[0];
x q[1];
rz(2.8802323) q[2];
sx q[2];
rz(-2.0730505) q[2];
sx q[2];
rz(2.4136191) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5708872) q[1];
sx q[1];
rz(-1.0859617) q[1];
sx q[1];
rz(0.2356694) q[1];
rz(0.44002779) q[3];
sx q[3];
rz(-1.3413359) q[3];
sx q[3];
rz(3.0408531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.51776) q[2];
sx q[2];
rz(-1.0972923) q[2];
sx q[2];
rz(0.1499873) q[2];
rz(-3.127626) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(2.6278031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(-0.81277043) q[0];
rz(1.4179519) q[1];
sx q[1];
rz(-1.5128472) q[1];
sx q[1];
rz(-2.0779804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3115754) q[0];
sx q[0];
rz(-0.63236134) q[0];
sx q[0];
rz(-0.57741369) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0399357) q[2];
sx q[2];
rz(-0.16451193) q[2];
sx q[2];
rz(-0.042334231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8001682) q[1];
sx q[1];
rz(-1.5593976) q[1];
sx q[1];
rz(-2.5428548) q[1];
rz(-pi) q[2];
x q[2];
rz(2.205414) q[3];
sx q[3];
rz(-2.8293316) q[3];
sx q[3];
rz(-1.043751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26760462) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.7248636) q[2];
rz(3.0806165) q[3];
sx q[3];
rz(-0.91106001) q[3];
sx q[3];
rz(-1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42162073) q[0];
sx q[0];
rz(-2.3201729) q[0];
sx q[0];
rz(0.11909568) q[0];
rz(-0.47856092) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(2.4518769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0473339) q[0];
sx q[0];
rz(-2.2660997) q[0];
sx q[0];
rz(-0.78710763) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1999667) q[2];
sx q[2];
rz(-1.1111819) q[2];
sx q[2];
rz(0.72061611) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9097773) q[1];
sx q[1];
rz(-1.6157774) q[1];
sx q[1];
rz(-0.86391967) q[1];
x q[2];
rz(-1.7411362) q[3];
sx q[3];
rz(-0.719845) q[3];
sx q[3];
rz(-0.094739044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0488284) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(1.0678585) q[2];
rz(-0.16718665) q[3];
sx q[3];
rz(-1.3719631) q[3];
sx q[3];
rz(-3.0021477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329413) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(2.2331878) q[0];
rz(-2.1482229) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(-0.87432528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8876111) q[0];
sx q[0];
rz(-1.7029666) q[0];
sx q[0];
rz(1.553276) q[0];
rz(1.3130593) q[2];
sx q[2];
rz(-1.589587) q[2];
sx q[2];
rz(-2.7692139) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0093551) q[1];
sx q[1];
rz(-1.7061846) q[1];
sx q[1];
rz(-3.0817249) q[1];
rz(-pi) q[2];
rz(2.1217974) q[3];
sx q[3];
rz(-1.2606818) q[3];
sx q[3];
rz(1.9581025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5270762) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(0.85421526) q[2];
rz(0.45520511) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(-0.2555041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47827569) q[0];
sx q[0];
rz(-1.217696) q[0];
sx q[0];
rz(1.664337) q[0];
rz(1.5669589) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(2.4050567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2701243) q[0];
sx q[0];
rz(-2.0379319) q[0];
sx q[0];
rz(0.26345912) q[0];
x q[1];
rz(-2.7608654) q[2];
sx q[2];
rz(-1.6672616) q[2];
sx q[2];
rz(-0.42489932) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7875576) q[1];
sx q[1];
rz(-1.8144467) q[1];
sx q[1];
rz(0.32275782) q[1];
rz(-pi) q[2];
rz(-0.035321354) q[3];
sx q[3];
rz(-0.97408726) q[3];
sx q[3];
rz(-0.91340706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7443202) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.5273904) q[2];
rz(-2.7496036) q[3];
sx q[3];
rz(-1.1992998) q[3];
sx q[3];
rz(0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024260661) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(-0.20183739) q[0];
rz(-0.2233389) q[1];
sx q[1];
rz(-2.6170862) q[1];
sx q[1];
rz(-2.7490659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46014547) q[0];
sx q[0];
rz(-2.0325615) q[0];
sx q[0];
rz(3.0236834) q[0];
rz(-pi) q[1];
rz(1.6286276) q[2];
sx q[2];
rz(-2.4041345) q[2];
sx q[2];
rz(-0.89227572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.10822693) q[1];
sx q[1];
rz(-1.997233) q[1];
sx q[1];
rz(2.1436177) q[1];
rz(-pi) q[2];
rz(1.2667381) q[3];
sx q[3];
rz(-2.0342329) q[3];
sx q[3];
rz(-0.26362905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4203804) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(0.58050275) q[2];
rz(0.37603363) q[3];
sx q[3];
rz(-2.1707363) q[3];
sx q[3];
rz(1.5413126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15039438) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(1.6509542) q[1];
sx q[1];
rz(-1.8883659) q[1];
sx q[1];
rz(-2.2813003) q[1];
rz(0.80581325) q[2];
sx q[2];
rz(-2.056682) q[2];
sx q[2];
rz(-0.32120612) q[2];
rz(-0.069546139) q[3];
sx q[3];
rz(-2.8290717) q[3];
sx q[3];
rz(1.1858218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
