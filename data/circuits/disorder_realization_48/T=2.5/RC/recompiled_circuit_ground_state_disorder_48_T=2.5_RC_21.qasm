OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1005062) q[0];
sx q[0];
rz(-2.9334928) q[0];
sx q[0];
rz(-2.7383374) q[0];
rz(-0.23303214) q[1];
sx q[1];
rz(-1.4401399) q[1];
sx q[1];
rz(-2.9177102) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96432811) q[0];
sx q[0];
rz(-1.3601662) q[0];
sx q[0];
rz(-0.62646477) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86647948) q[2];
sx q[2];
rz(-1.6687376) q[2];
sx q[2];
rz(2.3549454) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63131166) q[1];
sx q[1];
rz(-1.3461539) q[1];
sx q[1];
rz(-2.2366877) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86312189) q[3];
sx q[3];
rz(-0.89608367) q[3];
sx q[3];
rz(2.5503039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69748059) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(-3.1057788) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(2.9353976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8908454) q[0];
sx q[0];
rz(-2.146281) q[0];
sx q[0];
rz(-0.19749755) q[0];
rz(2.6644871) q[1];
sx q[1];
rz(-2.4767866) q[1];
sx q[1];
rz(0.8055996) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2458513) q[0];
sx q[0];
rz(-2.4516458) q[0];
sx q[0];
rz(0.2941546) q[0];
rz(-pi) q[1];
rz(-0.16519157) q[2];
sx q[2];
rz(-0.38536638) q[2];
sx q[2];
rz(-1.2549653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2195314) q[1];
sx q[1];
rz(-2.1967255) q[1];
sx q[1];
rz(-2.3971167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74975852) q[3];
sx q[3];
rz(-0.20515144) q[3];
sx q[3];
rz(1.6501901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8311367) q[2];
sx q[2];
rz(-1.4996935) q[2];
sx q[2];
rz(1.7384701) q[2];
rz(0.63203114) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(0.12701756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8359351) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(2.054731) q[0];
rz(0.19733363) q[1];
sx q[1];
rz(-1.6355762) q[1];
sx q[1];
rz(2.8089583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5088288) q[0];
sx q[0];
rz(-0.45464215) q[0];
sx q[0];
rz(-1.2484545) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.640075) q[2];
sx q[2];
rz(-2.4054619) q[2];
sx q[2];
rz(-0.87321216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3217993) q[1];
sx q[1];
rz(-1.7463356) q[1];
sx q[1];
rz(1.3513397) q[1];
x q[2];
rz(-2.4147291) q[3];
sx q[3];
rz(-1.3976025) q[3];
sx q[3];
rz(-1.7693335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7670333) q[2];
sx q[2];
rz(-1.550753) q[2];
sx q[2];
rz(1.6732875) q[2];
rz(3.1324006) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(2.3495242) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.511635) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(-1.4427503) q[0];
rz(1.0244145) q[1];
sx q[1];
rz(-1.2316278) q[1];
sx q[1];
rz(-0.82694298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073155135) q[0];
sx q[0];
rz(-0.19864635) q[0];
sx q[0];
rz(-1.9563) q[0];
rz(2.9873136) q[2];
sx q[2];
rz(-2.2529229) q[2];
sx q[2];
rz(2.9088504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0961541) q[1];
sx q[1];
rz(-1.7806345) q[1];
sx q[1];
rz(0.44053284) q[1];
rz(-pi) q[2];
rz(-0.8884646) q[3];
sx q[3];
rz(-1.4766221) q[3];
sx q[3];
rz(-1.7674131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3247165) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(1.8640222) q[2];
rz(-2.4534524) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-2.1791606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0012896) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-0.21743123) q[0];
rz(-0.28542074) q[1];
sx q[1];
rz(-0.60573429) q[1];
sx q[1];
rz(-0.4981471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.065154) q[0];
sx q[0];
rz(-0.74238837) q[0];
sx q[0];
rz(1.0774633) q[0];
rz(-0.26136036) q[2];
sx q[2];
rz(-1.0685421) q[2];
sx q[2];
rz(-2.4136191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5707055) q[1];
sx q[1];
rz(-1.0859617) q[1];
sx q[1];
rz(2.9059232) q[1];
rz(-pi) q[2];
rz(-0.44002779) q[3];
sx q[3];
rz(-1.3413359) q[3];
sx q[3];
rz(0.10073951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.51776) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(-0.1499873) q[2];
rz(3.127626) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(0.51378957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496395) q[0];
sx q[0];
rz(-0.54170251) q[0];
sx q[0];
rz(0.81277043) q[0];
rz(1.4179519) q[1];
sx q[1];
rz(-1.5128472) q[1];
sx q[1];
rz(-2.0779804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83001721) q[0];
sx q[0];
rz(-0.63236134) q[0];
sx q[0];
rz(-2.564179) q[0];
rz(-0.16367775) q[2];
sx q[2];
rz(-1.5874169) q[2];
sx q[2];
rz(-1.7134242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8955161) q[1];
sx q[1];
rz(-0.59883307) q[1];
sx q[1];
rz(3.1213698) q[1];
rz(-2.952488) q[3];
sx q[3];
rz(-1.8207887) q[3];
sx q[3];
rz(0.3853021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26760462) q[2];
sx q[2];
rz(-0.559811) q[2];
sx q[2];
rz(1.416729) q[2];
rz(-3.0806165) q[3];
sx q[3];
rz(-0.91106001) q[3];
sx q[3];
rz(1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.42162073) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-3.022497) q[0];
rz(0.47856092) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(0.68971577) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0473339) q[0];
sx q[0];
rz(-0.87549291) q[0];
sx q[0];
rz(-0.78710763) q[0];
rz(-2.5096171) q[2];
sx q[2];
rz(-2.5594829) q[2];
sx q[2];
rz(1.7012973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3773681) q[1];
sx q[1];
rz(-0.86478327) q[1];
sx q[1];
rz(-0.059125916) q[1];
x q[2];
rz(2.2834217) q[3];
sx q[3];
rz(-1.6827876) q[3];
sx q[3];
rz(-1.5369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0488284) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(-1.0678585) q[2];
rz(-0.16718665) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(-0.13944496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(-2.2331878) q[0];
rz(-2.1482229) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(-0.87432528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8876111) q[0];
sx q[0];
rz(-1.7029666) q[0];
sx q[0];
rz(-1.553276) q[0];
rz(-1.6443917) q[2];
sx q[2];
rz(-2.8831867) q[2];
sx q[2];
rz(-1.1272573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5501551) q[1];
sx q[1];
rz(-0.14796013) q[1];
sx q[1];
rz(1.9846538) q[1];
x q[2];
rz(-0.35975144) q[3];
sx q[3];
rz(-2.092741) q[3];
sx q[3];
rz(2.9396543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61451644) q[2];
sx q[2];
rz(-2.9824342) q[2];
sx q[2];
rz(-0.85421526) q[2];
rz(-0.45520511) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(-2.8860886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.663317) q[0];
sx q[0];
rz(-1.217696) q[0];
sx q[0];
rz(1.664337) q[0];
rz(-1.5669589) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(-2.4050567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5631396) q[0];
sx q[0];
rz(-1.805465) q[0];
sx q[0];
rz(-1.0893954) q[0];
rz(-pi) q[1];
rz(-2.8868448) q[2];
sx q[2];
rz(-2.7494135) q[2];
sx q[2];
rz(-2.2318411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29727259) q[1];
sx q[1];
rz(-1.2579009) q[1];
sx q[1];
rz(1.827153) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1677954) q[3];
sx q[3];
rz(-1.5415808) q[3];
sx q[3];
rz(0.67724281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3972724) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(-1.5273904) q[2];
rz(0.39198908) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(-0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024260661) q[0];
sx q[0];
rz(-1.6680822) q[0];
sx q[0];
rz(0.20183739) q[0];
rz(2.9182538) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(2.7490659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4215764) q[0];
sx q[0];
rz(-2.6660653) q[0];
sx q[0];
rz(-1.3386734) q[0];
rz(3.0891339) q[2];
sx q[2];
rz(-2.3067368) q[2];
sx q[2];
rz(-0.8142161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9397247) q[1];
sx q[1];
rz(-2.0868667) q[1];
sx q[1];
rz(0.49560541) q[1];
rz(-pi) q[2];
rz(2.6018312) q[3];
sx q[3];
rz(-0.54815147) q[3];
sx q[3];
rz(-0.87566151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4203804) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(-2.5610899) q[2];
rz(-2.765559) q[3];
sx q[3];
rz(-0.97085634) q[3];
sx q[3];
rz(-1.5413126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.15039438) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(1.4906384) q[1];
sx q[1];
rz(-1.2532267) q[1];
sx q[1];
rz(0.86029235) q[1];
rz(2.5096624) q[2];
sx q[2];
rz(-0.9117374) q[2];
sx q[2];
rz(0.82814817) q[2];
rz(1.5483472) q[3];
sx q[3];
rz(-1.2590564) q[3];
sx q[3];
rz(1.1127478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
