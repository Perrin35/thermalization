OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8695495) q[0];
sx q[0];
rz(-0.46177319) q[0];
sx q[0];
rz(0.053134993) q[0];
rz(2.7825836) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(-0.63149482) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8780958) q[1];
sx q[1];
rz(-0.43259183) q[1];
sx q[1];
rz(-2.2357975) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9840368) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(-0.36039263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(1.9155496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81235028) q[0];
sx q[0];
rz(-1.7698405) q[0];
sx q[0];
rz(0.45254405) q[0];
rz(-pi) q[1];
rz(-2.3833582) q[2];
sx q[2];
rz(-0.97823921) q[2];
sx q[2];
rz(2.638608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.023149816) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(-1.5688194) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5561043) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(-3.1327914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91988504) q[0];
sx q[0];
rz(-1.3516597) q[0];
sx q[0];
rz(2.1973781) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8573895) q[2];
sx q[2];
rz(-0.9409875) q[2];
sx q[2];
rz(2.0558002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2468977) q[1];
sx q[1];
rz(-2.9322335) q[1];
sx q[1];
rz(-0.73114242) q[1];
rz(-3.0089278) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-0.69561361) q[0];
sx q[0];
rz(0.22380933) q[0];
rz(0.085725351) q[2];
sx q[2];
rz(-1.013373) q[2];
sx q[2];
rz(2.2103708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.013853156) q[1];
sx q[1];
rz(-1.2989559) q[1];
sx q[1];
rz(-2.4747162) q[1];
rz(-pi) q[2];
rz(2.8193072) q[3];
sx q[3];
rz(-2.5282113) q[3];
sx q[3];
rz(-1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.5198583) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653523) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-0.68666896) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72563719) q[2];
sx q[2];
rz(-2.0113532) q[2];
sx q[2];
rz(0.63647905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4836854) q[1];
sx q[1];
rz(-0.92805082) q[1];
sx q[1];
rz(-2.583858) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(-2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.58364761) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(0.46498743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2705921) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(0.90765783) q[0];
rz(-pi) q[1];
rz(0.20435135) q[2];
sx q[2];
rz(-0.34826476) q[2];
sx q[2];
rz(-1.5262926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8556559) q[1];
sx q[1];
rz(-2.2054513) q[1];
sx q[1];
rz(-0.059307701) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3607849) q[3];
sx q[3];
rz(-2.2689153) q[3];
sx q[3];
rz(-0.39892808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(2.6965551) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(2.8745108) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(-3.1177915) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10931817) q[0];
sx q[0];
rz(-2.8335857) q[0];
sx q[0];
rz(-0.59662915) q[0];
rz(-pi) q[1];
rz(-2.6290226) q[2];
sx q[2];
rz(-1.1935496) q[2];
sx q[2];
rz(2.9104428) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1289039) q[1];
sx q[1];
rz(-1.745599) q[1];
sx q[1];
rz(0.054794475) q[1];
rz(-2.6050623) q[3];
sx q[3];
rz(-1.1324258) q[3];
sx q[3];
rz(-0.28552548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(2.8038213) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-0.24857323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08911207) q[0];
sx q[0];
rz(-0.57894527) q[0];
sx q[0];
rz(-0.46778932) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0806662) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(0.61444297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92330248) q[1];
sx q[1];
rz(-0.81070886) q[1];
sx q[1];
rz(0.87358012) q[1];
rz(-0.076898889) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.988525) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(-1.3185906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0193664) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(-3.0526572) q[0];
x q[1];
rz(-2.4852072) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(-0.93751794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1773771) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(-1.7978976) q[1];
x q[2];
rz(-1.1845469) q[3];
sx q[3];
rz(-2.2830314) q[3];
sx q[3];
rz(1.4854747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(0.88360751) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(-1.5426853) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0364089) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(1.4559792) q[0];
rz(-pi) q[1];
rz(-3.10896) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(2.0805217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1120158) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(-0.13346787) q[1];
rz(0.19461467) q[3];
sx q[3];
rz(-0.6932887) q[3];
sx q[3];
rz(2.4926536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-1.6115887) q[2];
sx q[2];
rz(-2.5695124) q[2];
sx q[2];
rz(3.1281501) q[2];
rz(-0.59109296) q[3];
sx q[3];
rz(-1.4553087) q[3];
sx q[3];
rz(1.6607264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];