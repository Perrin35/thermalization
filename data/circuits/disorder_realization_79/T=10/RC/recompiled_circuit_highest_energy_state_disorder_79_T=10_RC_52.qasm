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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(-0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(-1.5348943) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255768) q[0];
sx q[0];
rz(-1.1281779) q[0];
sx q[0];
rz(-0.26066633) q[0];
rz(-1.9594749) q[2];
sx q[2];
rz(-0.76180327) q[2];
sx q[2];
rz(-2.3488059) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.038329934) q[1];
sx q[1];
rz(-0.94417773) q[1];
sx q[1];
rz(0.4396529) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4356218) q[3];
sx q[3];
rz(-1.7723284) q[3];
sx q[3];
rz(0.38231787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(0.99754769) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(-1.0409748) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(1.8541699) q[0];
rz(2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.9226673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78271237) q[0];
sx q[0];
rz(-0.25280372) q[0];
sx q[0];
rz(0.88835277) q[0];
rz(1.7886247) q[2];
sx q[2];
rz(-2.2309897) q[2];
sx q[2];
rz(0.24581395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87172958) q[1];
sx q[1];
rz(-0.31530373) q[1];
sx q[1];
rz(1.0878776) q[1];
x q[2];
rz(-0.35393012) q[3];
sx q[3];
rz(-2.2316389) q[3];
sx q[3];
rz(-0.61001813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(0.62180579) q[2];
rz(-0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(0.59858876) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-0.52282202) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110958) q[0];
sx q[0];
rz(-1.8930264) q[0];
sx q[0];
rz(2.8186174) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36575138) q[2];
sx q[2];
rz(-2.7257082) q[2];
sx q[2];
rz(-1.6037841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8717958) q[1];
sx q[1];
rz(-1.7933473) q[1];
sx q[1];
rz(-2.6128164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9013635) q[3];
sx q[3];
rz(-1.2798314) q[3];
sx q[3];
rz(-1.3823079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(1.5175021) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(1.602518) q[0];
rz(-1.0333565) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(-1.296952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596969) q[0];
sx q[0];
rz(-2.2821626) q[0];
sx q[0];
rz(2.864896) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6455075) q[2];
sx q[2];
rz(-2.4444408) q[2];
sx q[2];
rz(-2.8139204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9819543) q[1];
sx q[1];
rz(-2.5667563) q[1];
sx q[1];
rz(-2.3153618) q[1];
rz(-pi) q[2];
rz(-0.4850895) q[3];
sx q[3];
rz(-2.5311433) q[3];
sx q[3];
rz(-1.4097253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(-2.0959334) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3359208) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(-0.51934284) q[0];
rz(-2.1995811) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(1.6090144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518884) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(2.4469482) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18073323) q[2];
sx q[2];
rz(-1.7109979) q[2];
sx q[2];
rz(3.0716346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0253801) q[1];
sx q[1];
rz(-1.0702225) q[1];
sx q[1];
rz(-3.1356372) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4306253) q[3];
sx q[3];
rz(-2.9971854) q[3];
sx q[3];
rz(-2.292423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(-2.8572148) q[2];
rz(2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(-2.8936581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(0.64055842) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-2.9277149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.36422) q[0];
sx q[0];
rz(-1.7988482) q[0];
sx q[0];
rz(1.26012) q[0];
x q[1];
rz(1.1919695) q[2];
sx q[2];
rz(-2.2791499) q[2];
sx q[2];
rz(1.5713991) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8941514) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(3.0559866) q[1];
rz(-2.0908666) q[3];
sx q[3];
rz(-1.5087391) q[3];
sx q[3];
rz(-0.37953916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(2.9638929) q[2];
rz(-1.3972345) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(-2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099667065) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(-1.1055111) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.6800605) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34237) q[0];
sx q[0];
rz(-0.91081753) q[0];
sx q[0];
rz(2.2837385) q[0];
rz(-pi) q[1];
rz(-2.2093532) q[2];
sx q[2];
rz(-1.1503714) q[2];
sx q[2];
rz(2.1115542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.098432861) q[1];
sx q[1];
rz(-2.5665356) q[1];
sx q[1];
rz(0.21456031) q[1];
rz(1.8959778) q[3];
sx q[3];
rz(-1.0187314) q[3];
sx q[3];
rz(-1.0528885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(2.936787) q[2];
rz(-1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.4749984) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(-1.9974476) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0856253) q[0];
sx q[0];
rz(-0.48729839) q[0];
sx q[0];
rz(3.0874599) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2828243) q[2];
sx q[2];
rz(-0.33703732) q[2];
sx q[2];
rz(-2.2971643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.543964) q[1];
sx q[1];
rz(-1.9457327) q[1];
sx q[1];
rz(-2.4191394) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89544483) q[3];
sx q[3];
rz(-2.709528) q[3];
sx q[3];
rz(2.2275782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7555776) q[2];
sx q[2];
rz(-2.3365946) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-2.2224865) q[3];
sx q[3];
rz(-0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(1.8852604) q[0];
rz(-0.095257692) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(0.5307861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7266709) q[0];
sx q[0];
rz(-1.4106531) q[0];
sx q[0];
rz(-1.8216351) q[0];
x q[1];
rz(-1.5321391) q[2];
sx q[2];
rz(-1.8886023) q[2];
sx q[2];
rz(0.048887756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8618882) q[1];
sx q[1];
rz(-1.8100396) q[1];
sx q[1];
rz(1.4390535) q[1];
rz(-pi) q[2];
rz(2.9050499) q[3];
sx q[3];
rz(-2.6841087) q[3];
sx q[3];
rz(1.2381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(-2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(2.7428108) q[0];
rz(-1.3840236) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(-0.049364518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684262) q[0];
sx q[0];
rz(-0.12463649) q[0];
sx q[0];
rz(-1.0602555) q[0];
rz(-1.1902383) q[2];
sx q[2];
rz(-2.8597288) q[2];
sx q[2];
rz(-0.50279248) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70061848) q[1];
sx q[1];
rz(-1.6153533) q[1];
sx q[1];
rz(-0.99127165) q[1];
rz(-pi) q[2];
rz(-1.7743763) q[3];
sx q[3];
rz(-0.94881159) q[3];
sx q[3];
rz(-0.79727117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-1.0646959) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(-0.76416107) q[2];
sx q[2];
rz(-2.4780826) q[2];
sx q[2];
rz(1.6586951) q[2];
rz(-1.0743027) q[3];
sx q[3];
rz(-2.5357694) q[3];
sx q[3];
rz(2.3800935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
