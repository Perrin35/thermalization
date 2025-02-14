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
rz(-0.66145575) q[0];
sx q[0];
rz(2.8602726) q[0];
sx q[0];
rz(8.2814132) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(-1.3209359) q[1];
sx q[1];
rz(0.015425711) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93958873) q[0];
sx q[0];
rz(-2.1103854) q[0];
sx q[0];
rz(2.5651776) q[0];
x q[1];
rz(0.5792497) q[2];
sx q[2];
rz(-2.0655895) q[2];
sx q[2];
rz(0.65828568) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4651067) q[1];
sx q[1];
rz(-0.78395212) q[1];
sx q[1];
rz(2.1539861) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3728834) q[3];
sx q[3];
rz(-2.1002401) q[3];
sx q[3];
rz(2.6967665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.780484) q[2];
sx q[2];
rz(-0.0094272308) q[2];
sx q[2];
rz(1.1040322) q[2];
rz(-0.55936724) q[3];
sx q[3];
rz(-2.1909824) q[3];
sx q[3];
rz(2.2573788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169432) q[0];
sx q[0];
rz(-2.556612) q[0];
sx q[0];
rz(-0.90644932) q[0];
rz(-0.34140423) q[1];
sx q[1];
rz(-0.87541348) q[1];
sx q[1];
rz(-0.34814775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6432836) q[0];
sx q[0];
rz(-1.3237244) q[0];
sx q[0];
rz(2.8227795) q[0];
x q[1];
rz(1.6115336) q[2];
sx q[2];
rz(-1.5447283) q[2];
sx q[2];
rz(3.0422826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9056919) q[1];
sx q[1];
rz(-2.5254619) q[1];
sx q[1];
rz(-1.4409164) q[1];
x q[2];
rz(0.82795469) q[3];
sx q[3];
rz(-1.7418234) q[3];
sx q[3];
rz(-3.1307182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1657095) q[2];
sx q[2];
rz(-0.3419064) q[2];
sx q[2];
rz(0.085414097) q[2];
rz(-0.092770569) q[3];
sx q[3];
rz(-2.2723891) q[3];
sx q[3];
rz(-2.2658277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669325) q[0];
sx q[0];
rz(-0.36151883) q[0];
sx q[0];
rz(0.34459484) q[0];
rz(1.5498281) q[1];
sx q[1];
rz(-1.7235618) q[1];
sx q[1];
rz(1.6835015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8710468) q[0];
sx q[0];
rz(-0.71840175) q[0];
sx q[0];
rz(-3.0860391) q[0];
rz(0.05942101) q[2];
sx q[2];
rz(-2.0869617) q[2];
sx q[2];
rz(-2.8490861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5815253) q[1];
sx q[1];
rz(-1.4017644) q[1];
sx q[1];
rz(1.8174816) q[1];
x q[2];
rz(-0.97403861) q[3];
sx q[3];
rz(-1.1263873) q[3];
sx q[3];
rz(0.35154464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6452667) q[2];
sx q[2];
rz(-2.2023449) q[2];
sx q[2];
rz(-1.1305031) q[2];
rz(-0.94275236) q[3];
sx q[3];
rz(-2.0944984) q[3];
sx q[3];
rz(-1.2070791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4261674) q[0];
sx q[0];
rz(-1.796145) q[0];
sx q[0];
rz(1.4581534) q[0];
rz(2.093105) q[1];
sx q[1];
rz(-1.6494992) q[1];
sx q[1];
rz(2.6715211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78280503) q[0];
sx q[0];
rz(-2.067152) q[0];
sx q[0];
rz(1.8904782) q[0];
rz(-1.1328892) q[2];
sx q[2];
rz(-0.57269579) q[2];
sx q[2];
rz(-1.1066135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25950228) q[1];
sx q[1];
rz(-0.90356245) q[1];
sx q[1];
rz(-1.3253071) q[1];
rz(-pi) q[2];
rz(-0.87155452) q[3];
sx q[3];
rz(-2.1208753) q[3];
sx q[3];
rz(-0.14000237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5992392) q[2];
sx q[2];
rz(-2.3848644) q[2];
sx q[2];
rz(-0.92764628) q[2];
rz(-1.2841691) q[3];
sx q[3];
rz(-1.7312867) q[3];
sx q[3];
rz(-0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014378431) q[0];
sx q[0];
rz(-0.40507409) q[0];
sx q[0];
rz(2.6601484) q[0];
rz(-2.1638347) q[1];
sx q[1];
rz(-0.6441741) q[1];
sx q[1];
rz(2.0668623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.510493) q[0];
sx q[0];
rz(-1.298615) q[0];
sx q[0];
rz(-1.136632) q[0];
rz(1.2891026) q[2];
sx q[2];
rz(-1.8605202) q[2];
sx q[2];
rz(-0.85838028) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2664484) q[1];
sx q[1];
rz(-1.9423264) q[1];
sx q[1];
rz(-2.8846011) q[1];
rz(-pi) q[2];
rz(2.7513192) q[3];
sx q[3];
rz(-1.8650132) q[3];
sx q[3];
rz(-0.38021429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0398728) q[2];
sx q[2];
rz(-1.1168672) q[2];
sx q[2];
rz(-2.5679585) q[2];
rz(1.6576069) q[3];
sx q[3];
rz(-2.0935757) q[3];
sx q[3];
rz(-0.60605961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023733519) q[0];
sx q[0];
rz(-2.8890299) q[0];
sx q[0];
rz(-2.8299487) q[0];
rz(0.32870865) q[1];
sx q[1];
rz(-1.3361822) q[1];
sx q[1];
rz(2.4005344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8097654) q[0];
sx q[0];
rz(-2.1719195) q[0];
sx q[0];
rz(-2.580216) q[0];
rz(1.8297878) q[2];
sx q[2];
rz(-0.2364859) q[2];
sx q[2];
rz(-1.8571929) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.44081008) q[1];
sx q[1];
rz(-0.68432552) q[1];
sx q[1];
rz(-1.6709063) q[1];
rz(-pi) q[2];
rz(0.00087329207) q[3];
sx q[3];
rz(-0.93993087) q[3];
sx q[3];
rz(1.6959977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3662423) q[2];
sx q[2];
rz(-2.5516208) q[2];
sx q[2];
rz(0.74014202) q[2];
rz(-0.38672334) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(2.6630785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1665523) q[0];
sx q[0];
rz(-0.98169011) q[0];
sx q[0];
rz(-2.8983086) q[0];
rz(0.89726204) q[1];
sx q[1];
rz(-1.5586531) q[1];
sx q[1];
rz(-2.3975587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720111) q[0];
sx q[0];
rz(-1.9599304) q[0];
sx q[0];
rz(1.6604108) q[0];
rz(-pi) q[1];
rz(2.7626708) q[2];
sx q[2];
rz(-2.1728467) q[2];
sx q[2];
rz(-0.10720358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0571088) q[1];
sx q[1];
rz(-1.9388102) q[1];
sx q[1];
rz(-2.1294566) q[1];
rz(-0.18455581) q[3];
sx q[3];
rz(-1.0519487) q[3];
sx q[3];
rz(-2.542582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1870785) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(-3.1336866) q[2];
rz(-1.6451969) q[3];
sx q[3];
rz(-3.0683066) q[3];
sx q[3];
rz(2.6325398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1169443) q[0];
sx q[0];
rz(-3.109303) q[0];
sx q[0];
rz(3.0049128) q[0];
rz(-1.7928803) q[1];
sx q[1];
rz(-1.8486479) q[1];
sx q[1];
rz(2.6332556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0011372066) q[0];
sx q[0];
rz(-1.1759725) q[0];
sx q[0];
rz(0.14060256) q[0];
x q[1];
rz(-2.2466564) q[2];
sx q[2];
rz(-1.8064665) q[2];
sx q[2];
rz(-1.6431944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4813684) q[1];
sx q[1];
rz(-1.9268039) q[1];
sx q[1];
rz(-2.2769209) q[1];
rz(-pi) q[2];
rz(-1.6898722) q[3];
sx q[3];
rz(-1.4919466) q[3];
sx q[3];
rz(-0.54311968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4623744) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(0.057961658) q[2];
rz(-2.1884749) q[3];
sx q[3];
rz(-2.0942196) q[3];
sx q[3];
rz(-0.63547772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304831) q[0];
sx q[0];
rz(-1.4169175) q[0];
sx q[0];
rz(-0.86607754) q[0];
rz(1.3653612) q[1];
sx q[1];
rz(-1.5581286) q[1];
sx q[1];
rz(-2.3376215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3063076) q[0];
sx q[0];
rz(-0.87125766) q[0];
sx q[0];
rz(0.37658306) q[0];
rz(-pi) q[1];
rz(2.3525904) q[2];
sx q[2];
rz(-1.2478634) q[2];
sx q[2];
rz(-0.6168405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2197847) q[1];
sx q[1];
rz(-0.91118813) q[1];
sx q[1];
rz(-0.54306026) q[1];
rz(-pi) q[2];
rz(1.6220408) q[3];
sx q[3];
rz(-0.4040904) q[3];
sx q[3];
rz(1.3431637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0552401) q[2];
sx q[2];
rz(-0.1487727) q[2];
sx q[2];
rz(-1.0521592) q[2];
rz(2.2822028) q[3];
sx q[3];
rz(-2.342577) q[3];
sx q[3];
rz(-0.94256443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6732442) q[0];
sx q[0];
rz(-0.66272658) q[0];
sx q[0];
rz(0.41917875) q[0];
rz(0.016050922) q[1];
sx q[1];
rz(-1.5546067) q[1];
sx q[1];
rz(0.12414653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57174411) q[0];
sx q[0];
rz(-2.1198065) q[0];
sx q[0];
rz(0.17904539) q[0];
x q[1];
rz(-2.2474278) q[2];
sx q[2];
rz(-0.32197201) q[2];
sx q[2];
rz(-1.3077259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.032822306) q[1];
sx q[1];
rz(-1.747393) q[1];
sx q[1];
rz(-2.887421) q[1];
rz(-1.8087555) q[3];
sx q[3];
rz(-0.66413022) q[3];
sx q[3];
rz(-2.7414764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19196709) q[2];
sx q[2];
rz(-2.4032335) q[2];
sx q[2];
rz(2.9730566) q[2];
rz(1.656823) q[3];
sx q[3];
rz(-2.6717581) q[3];
sx q[3];
rz(-2.7114939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94854245) q[0];
sx q[0];
rz(-0.47946231) q[0];
sx q[0];
rz(-0.23647501) q[0];
rz(-0.82462689) q[1];
sx q[1];
rz(-1.5390479) q[1];
sx q[1];
rz(1.8631757) q[1];
rz(0.77565754) q[2];
sx q[2];
rz(-1.1393329) q[2];
sx q[2];
rz(2.3928497) q[2];
rz(-1.0239368) q[3];
sx q[3];
rz(-1.5187859) q[3];
sx q[3];
rz(-2.8092629) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
