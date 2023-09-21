OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(3.9290805) q[1];
sx q[1];
rz(10.706283) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7377388) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(-1.2234729) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69252695) q[2];
sx q[2];
rz(-1.5934957) q[2];
sx q[2];
rz(1.2251309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(-1.185226) q[1];
rz(-pi) q[2];
rz(2.0251861) q[3];
sx q[3];
rz(-0.14305275) q[3];
sx q[3];
rz(-1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.981502) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(-1.7574969) q[0];
rz(-pi) q[1];
rz(2.312617) q[2];
sx q[2];
rz(-1.0114705) q[2];
sx q[2];
rz(-1.0499357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21178791) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(2.2499229) q[1];
rz(-3.00499) q[3];
sx q[3];
rz(-1.6009637) q[3];
sx q[3];
rz(0.68340106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.9525607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3184568) q[0];
sx q[0];
rz(-1.4134777) q[0];
sx q[0];
rz(2.8419028) q[0];
rz(-1.5306773) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(-0.2085533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6310198) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(2.2651947) q[1];
rz(-pi) q[2];
rz(0.69840188) q[3];
sx q[3];
rz(-1.9365371) q[3];
sx q[3];
rz(-1.48667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90551463) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.8006515) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.9569424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957164) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(-1.7488259) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2934154) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.4959469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2733172) q[1];
sx q[1];
rz(-1.6874716) q[1];
sx q[1];
rz(2.1627746) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42479997) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(-2.6257023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(-2.284164) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(2.4127814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279014) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(-1.6708899) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9786733) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(-1.7999072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5908969) q[1];
sx q[1];
rz(-1.3662852) q[1];
sx q[1];
rz(1.2414316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2767931) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52299243) q[0];
sx q[0];
rz(-0.97609659) q[0];
sx q[0];
rz(2.3408875) q[0];
rz(-3.0350787) q[2];
sx q[2];
rz(-1.4918367) q[2];
sx q[2];
rz(0.67895141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8710528) q[1];
sx q[1];
rz(-0.45696837) q[1];
sx q[1];
rz(0.096710042) q[1];
rz(-pi) q[2];
rz(-0.75943767) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1202576) q[0];
sx q[0];
rz(-1.406167) q[0];
sx q[0];
rz(0.067532587) q[0];
rz(0.61755007) q[2];
sx q[2];
rz(-2.6744161) q[2];
sx q[2];
rz(0.10305603) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9387503) q[1];
sx q[1];
rz(-0.65333594) q[1];
sx q[1];
rz(-2.5696239) q[1];
rz(-pi) q[2];
rz(1.6129937) q[3];
sx q[3];
rz(-1.183126) q[3];
sx q[3];
rz(-1.8488415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(-1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(-0.45773488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477752) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(1.093822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55191314) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(-1.0362831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(1.5921028) q[1];
rz(-1.0192972) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(-2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.2618077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5926873) q[0];
sx q[0];
rz(-1.3015916) q[0];
sx q[0];
rz(-0.44434987) q[0];
x q[1];
rz(-1.3627522) q[2];
sx q[2];
rz(-2.2135452) q[2];
sx q[2];
rz(-0.3026697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-0.8926549) q[1];
sx q[1];
rz(-1.2441586) q[1];
rz(-pi) q[2];
rz(-1.3654361) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(-0.53232771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(-0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.2465729) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(1.0704401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265855) q[0];
sx q[0];
rz(-0.75782776) q[0];
sx q[0];
rz(-1.3269781) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.038254914) q[2];
sx q[2];
rz(-2.3352211) q[2];
sx q[2];
rz(0.56626608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4813303) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(0.34923133) q[1];
rz(-pi) q[2];
rz(-0.6108547) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(-0.50869298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76876172) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.6147511) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.6282745) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-3.1038531) q[2];
sx q[2];
rz(-0.83913091) q[2];
sx q[2];
rz(-0.057374949) q[2];
rz(2.3970849) q[3];
sx q[3];
rz(-0.67529882) q[3];
sx q[3];
rz(-2.8532367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];