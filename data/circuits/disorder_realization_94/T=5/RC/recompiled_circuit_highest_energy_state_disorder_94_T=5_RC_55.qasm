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
rz(0.87297451) q[0];
sx q[0];
rz(-0.31384808) q[0];
sx q[0];
rz(2.0883972) q[0];
rz(1.6246417) q[1];
sx q[1];
rz(-0.2927953) q[1];
sx q[1];
rz(0.93661469) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1952269) q[0];
sx q[0];
rz(-1.5527301) q[0];
sx q[0];
rz(1.0626331) q[0];
rz(2.1457315) q[2];
sx q[2];
rz(-2.323192) q[2];
sx q[2];
rz(0.40413293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17736152) q[1];
sx q[1];
rz(-1.8500016) q[1];
sx q[1];
rz(2.5482168) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39150146) q[3];
sx q[3];
rz(-1.3927407) q[3];
sx q[3];
rz(2.4181929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34095731) q[2];
sx q[2];
rz(-2.1850047) q[2];
sx q[2];
rz(1.0613649) q[2];
rz(-2.495885) q[3];
sx q[3];
rz(-2.784745) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(2.6820768) q[0];
rz(-1.3872321) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(2.0769108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6613805) q[0];
sx q[0];
rz(-0.35786942) q[0];
sx q[0];
rz(2.6741739) q[0];
rz(-pi) q[1];
x q[1];
rz(0.063195066) q[2];
sx q[2];
rz(-0.26775751) q[2];
sx q[2];
rz(-1.4227941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.766651) q[1];
sx q[1];
rz(-1.8901575) q[1];
sx q[1];
rz(1.6691895) q[1];
x q[2];
rz(0.27402409) q[3];
sx q[3];
rz(-1.8962057) q[3];
sx q[3];
rz(-0.39135483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3620944) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(2.8703459) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(2.0103683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75854492) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(2.7591822) q[0];
rz(2.8657939) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(2.3013505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58038515) q[0];
sx q[0];
rz(-1.9044283) q[0];
sx q[0];
rz(1.1949657) q[0];
rz(-pi) q[1];
rz(-0.20277299) q[2];
sx q[2];
rz(-1.2744858) q[2];
sx q[2];
rz(-1.2057829) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.099681252) q[1];
sx q[1];
rz(-0.96158342) q[1];
sx q[1];
rz(0.64263338) q[1];
x q[2];
rz(2.6557585) q[3];
sx q[3];
rz(-1.6260901) q[3];
sx q[3];
rz(3.0818617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8648839) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(-3.0806105) q[2];
rz(-1.3649155) q[3];
sx q[3];
rz(-0.18326062) q[3];
sx q[3];
rz(2.0257559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10114305) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-2.2818991) q[0];
rz(-2.1172093) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(-0.76622564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9767155) q[0];
sx q[0];
rz(-0.055792965) q[0];
sx q[0];
rz(2.4528422) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0417074) q[2];
sx q[2];
rz(-1.7211564) q[2];
sx q[2];
rz(0.1783812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1376201) q[1];
sx q[1];
rz(-1.3070053) q[1];
sx q[1];
rz(-0.54171087) q[1];
x q[2];
rz(-2.2732417) q[3];
sx q[3];
rz(-1.5159226) q[3];
sx q[3];
rz(-0.026264852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.031294558) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(0.026570126) q[2];
rz(2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37102315) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(-0.41719607) q[0];
rz(-2.5542651) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(2.3951098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74325753) q[0];
sx q[0];
rz(-1.8881838) q[0];
sx q[0];
rz(-0.36821915) q[0];
rz(-pi) q[1];
rz(2.339614) q[2];
sx q[2];
rz(-0.97014135) q[2];
sx q[2];
rz(-0.13896143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44297781) q[1];
sx q[1];
rz(-1.9483951) q[1];
sx q[1];
rz(2.5966132) q[1];
rz(-pi) q[2];
rz(-2.6559791) q[3];
sx q[3];
rz(-1.1385185) q[3];
sx q[3];
rz(1.2114482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31009659) q[2];
sx q[2];
rz(-2.8079171) q[2];
sx q[2];
rz(-0.92158544) q[2];
rz(1.0725675) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(-2.5785562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614586) q[0];
sx q[0];
rz(-2.7549094) q[0];
sx q[0];
rz(2.4899546) q[0];
rz(0.98546511) q[1];
sx q[1];
rz(-2.5990504) q[1];
sx q[1];
rz(0.14269565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7008377) q[0];
sx q[0];
rz(-1.6336226) q[0];
sx q[0];
rz(1.8453034) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8593117) q[2];
sx q[2];
rz(-2.9476894) q[2];
sx q[2];
rz(2.1197723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7276926) q[1];
sx q[1];
rz(-1.6928612) q[1];
sx q[1];
rz(2.6660566) q[1];
rz(-pi) q[2];
rz(-1.296259) q[3];
sx q[3];
rz(-2.6724788) q[3];
sx q[3];
rz(2.9363471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(2.209668) q[2];
rz(-0.43632397) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(1.1599734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.49509224) q[0];
sx q[0];
rz(-0.91630542) q[0];
sx q[0];
rz(-1.5402933) q[0];
rz(-2.1401999) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(-0.95759773) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0315899) q[0];
sx q[0];
rz(-0.45781198) q[0];
sx q[0];
rz(0.97719595) q[0];
x q[1];
rz(0.25156542) q[2];
sx q[2];
rz(-1.4770763) q[2];
sx q[2];
rz(-2.9523015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55064978) q[1];
sx q[1];
rz(-0.47524449) q[1];
sx q[1];
rz(-1.5859423) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5374588) q[3];
sx q[3];
rz(-1.0126202) q[3];
sx q[3];
rz(-0.5287002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6340948) q[2];
sx q[2];
rz(-1.2853421) q[2];
sx q[2];
rz(1.1320587) q[2];
rz(0.13103983) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(-1.8989782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2496654) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(3.0666572) q[0];
rz(0.19649188) q[1];
sx q[1];
rz(-0.34759977) q[1];
sx q[1];
rz(2.4620655) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.644548) q[0];
sx q[0];
rz(-1.0168075) q[0];
sx q[0];
rz(-1.8084256) q[0];
rz(-pi) q[1];
rz(-1.9747693) q[2];
sx q[2];
rz(-1.543664) q[2];
sx q[2];
rz(0.89744842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7767002) q[1];
sx q[1];
rz(-0.51799315) q[1];
sx q[1];
rz(-0.97229506) q[1];
rz(-pi) q[2];
rz(-1.5854168) q[3];
sx q[3];
rz(-1.059766) q[3];
sx q[3];
rz(2.8255812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3139265) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(-2.0619681) q[2];
rz(-0.45352724) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(2.1760904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249348) q[0];
sx q[0];
rz(-0.17913945) q[0];
sx q[0];
rz(-3.0058885) q[0];
rz(2.9756359) q[1];
sx q[1];
rz(-0.9981007) q[1];
sx q[1];
rz(-0.96806324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57763273) q[0];
sx q[0];
rz(-1.5819307) q[0];
sx q[0];
rz(2.0861113) q[0];
rz(-1.8161681) q[2];
sx q[2];
rz(-2.5858063) q[2];
sx q[2];
rz(-1.1664558) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8598946) q[1];
sx q[1];
rz(-1.7813695) q[1];
sx q[1];
rz(0.50986664) q[1];
x q[2];
rz(-0.37517199) q[3];
sx q[3];
rz(-1.811677) q[3];
sx q[3];
rz(-1.8011805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(-0.53259069) q[2];
rz(0.79832625) q[3];
sx q[3];
rz(-2.7440378) q[3];
sx q[3];
rz(-3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.3548729) q[0];
sx q[0];
rz(-0.69364554) q[0];
sx q[0];
rz(-0.68674809) q[0];
rz(1.2314388) q[1];
sx q[1];
rz(-1.8309007) q[1];
sx q[1];
rz(-0.062006921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69294482) q[0];
sx q[0];
rz(-0.8668859) q[0];
sx q[0];
rz(1.4848029) q[0];
rz(-pi) q[1];
rz(-0.36293928) q[2];
sx q[2];
rz(-2.937995) q[2];
sx q[2];
rz(-2.3516977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7534291) q[1];
sx q[1];
rz(-1.7723188) q[1];
sx q[1];
rz(2.9934197) q[1];
x q[2];
rz(1.5560772) q[3];
sx q[3];
rz(-1.2728661) q[3];
sx q[3];
rz(-3.0760279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-2.176602) q[2];
sx q[2];
rz(0.17267257) q[2];
rz(0.73575819) q[3];
sx q[3];
rz(-0.57307214) q[3];
sx q[3];
rz(0.47634038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19685766) q[0];
sx q[0];
rz(-1.6898962) q[0];
sx q[0];
rz(1.663399) q[0];
rz(-0.860515) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(2.4233107) q[2];
sx q[2];
rz(-0.16213308) q[2];
sx q[2];
rz(-2.8181974) q[2];
rz(0.92100704) q[3];
sx q[3];
rz(-2.3719127) q[3];
sx q[3];
rz(1.0815892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
