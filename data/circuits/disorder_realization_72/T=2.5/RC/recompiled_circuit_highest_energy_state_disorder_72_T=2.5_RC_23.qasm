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
rz(0.82062757) q[0];
sx q[0];
rz(2.4790915) q[0];
sx q[0];
rz(9.6883246) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(4.4237408) q[1];
sx q[1];
rz(10.10034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1206664) q[0];
sx q[0];
rz(-1.5772737) q[0];
sx q[0];
rz(1.2351888) q[0];
x q[1];
rz(-2.3604806) q[2];
sx q[2];
rz(-0.75133577) q[2];
sx q[2];
rz(1.2349811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88982595) q[1];
sx q[1];
rz(-1.7210362) q[1];
sx q[1];
rz(-1.9637693) q[1];
rz(-pi) q[2];
rz(-1.0682929) q[3];
sx q[3];
rz(-1.2811105) q[3];
sx q[3];
rz(2.009249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88966173) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(2.7737854) q[2];
rz(0.31428549) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(-2.9634326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3322068) q[0];
sx q[0];
rz(-3.0520913) q[0];
sx q[0];
rz(1.7820763) q[0];
rz(1.8571732) q[1];
sx q[1];
rz(-1.1651243) q[1];
sx q[1];
rz(-3.0899835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92049128) q[0];
sx q[0];
rz(-1.8360385) q[0];
sx q[0];
rz(-2.9539501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.935595) q[2];
sx q[2];
rz(-1.9340746) q[2];
sx q[2];
rz(1.3000803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9533128) q[1];
sx q[1];
rz(-1.1968391) q[1];
sx q[1];
rz(3.1409821) q[1];
rz(2.4447024) q[3];
sx q[3];
rz(-2.5514388) q[3];
sx q[3];
rz(-1.3499638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0525558) q[2];
sx q[2];
rz(-2.8812228) q[2];
sx q[2];
rz(-0.54711771) q[2];
rz(-1.268092) q[3];
sx q[3];
rz(-1.781446) q[3];
sx q[3];
rz(0.33856302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7608305) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(-1.8007675) q[0];
rz(-0.85397589) q[1];
sx q[1];
rz(-1.2869336) q[1];
sx q[1];
rz(-2.8391848) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6198928) q[0];
sx q[0];
rz(-0.95657737) q[0];
sx q[0];
rz(0.34003652) q[0];
x q[1];
rz(0.45291169) q[2];
sx q[2];
rz(-0.49395091) q[2];
sx q[2];
rz(0.043722186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21061347) q[1];
sx q[1];
rz(-2.180778) q[1];
sx q[1];
rz(0.90175866) q[1];
rz(-pi) q[2];
rz(-1.0875999) q[3];
sx q[3];
rz(-0.54387605) q[3];
sx q[3];
rz(-2.9915006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93620187) q[2];
sx q[2];
rz(-1.256029) q[2];
sx q[2];
rz(2.5993247) q[2];
rz(-0.99644709) q[3];
sx q[3];
rz(-0.22765972) q[3];
sx q[3];
rz(-3.0136133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19313136) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(1.8782072) q[0];
rz(0.43426934) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(-1.451452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5562486) q[0];
sx q[0];
rz(-2.5894269) q[0];
sx q[0];
rz(2.564173) q[0];
x q[1];
rz(-0.66777467) q[2];
sx q[2];
rz(-2.479907) q[2];
sx q[2];
rz(2.0354605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60148224) q[1];
sx q[1];
rz(-2.9992691) q[1];
sx q[1];
rz(-0.92592923) q[1];
x q[2];
rz(-2.1918815) q[3];
sx q[3];
rz(-1.1009645) q[3];
sx q[3];
rz(2.7787152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0749391) q[2];
sx q[2];
rz(-0.94330072) q[2];
sx q[2];
rz(-2.2554876) q[2];
rz(3.1352299) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(-1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56115443) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(1.2503257) q[0];
rz(0.23288947) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(0.09045352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0165774) q[0];
sx q[0];
rz(-2.8629486) q[0];
sx q[0];
rz(-1.11193) q[0];
rz(-pi) q[1];
rz(-2.6695741) q[2];
sx q[2];
rz(-1.3296122) q[2];
sx q[2];
rz(0.84005737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.131784) q[1];
sx q[1];
rz(-2.4854641) q[1];
sx q[1];
rz(3.0387278) q[1];
rz(1.2430322) q[3];
sx q[3];
rz(-1.48945) q[3];
sx q[3];
rz(1.7880754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6470486) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(1.9336644) q[2];
rz(2.2203994) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-0.45879656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0334466) q[0];
sx q[0];
rz(-1.4220081) q[0];
sx q[0];
rz(2.6500927) q[0];
rz(-2.4132701) q[1];
sx q[1];
rz(-1.1110543) q[1];
sx q[1];
rz(-0.19457766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3964243) q[0];
sx q[0];
rz(-2.130742) q[0];
sx q[0];
rz(-0.26225175) q[0];
x q[1];
rz(0.75677432) q[2];
sx q[2];
rz(-0.96454731) q[2];
sx q[2];
rz(-0.9204677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3272596) q[1];
sx q[1];
rz(-0.93686283) q[1];
sx q[1];
rz(2.4933706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.803502) q[3];
sx q[3];
rz(-0.92706087) q[3];
sx q[3];
rz(2.2542867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3538094) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(0.0018399012) q[2];
rz(-1.4932102) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(-0.28436896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(-1.0650241) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-0.5717259) q[1];
sx q[1];
rz(0.46331847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36998591) q[0];
sx q[0];
rz(-2.5328203) q[0];
sx q[0];
rz(1.2941525) q[0];
x q[1];
rz(-0.58775522) q[2];
sx q[2];
rz(-2.3238733) q[2];
sx q[2];
rz(0.14434926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3101497) q[1];
sx q[1];
rz(-2.0528704) q[1];
sx q[1];
rz(-1.4170013) q[1];
x q[2];
rz(-1.4557119) q[3];
sx q[3];
rz(-1.5466549) q[3];
sx q[3];
rz(2.3379282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.00039438417) q[2];
sx q[2];
rz(-1.1606777) q[2];
sx q[2];
rz(-2.311643) q[2];
rz(-2.0937894) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-0.1786264) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33775) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(-0.75482279) q[0];
rz(-2.7524475) q[1];
sx q[1];
rz(-1.4578338) q[1];
sx q[1];
rz(-2.7438502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5311814) q[0];
sx q[0];
rz(-1.2875778) q[0];
sx q[0];
rz(-0.2679268) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7806578) q[2];
sx q[2];
rz(-0.58141469) q[2];
sx q[2];
rz(1.5370528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3263742) q[1];
sx q[1];
rz(-2.7644016) q[1];
sx q[1];
rz(2.3227824) q[1];
rz(-0.74716178) q[3];
sx q[3];
rz(-2.0748169) q[3];
sx q[3];
rz(0.51591831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8336746) q[2];
sx q[2];
rz(-2.8454744) q[2];
sx q[2];
rz(2.8274242) q[2];
rz(-1.8999772) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(1.2392905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32387963) q[0];
sx q[0];
rz(-2.4535024) q[0];
sx q[0];
rz(-0.24539465) q[0];
rz(-1.2303526) q[1];
sx q[1];
rz(-3.0407258) q[1];
sx q[1];
rz(0.51602498) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9629899) q[0];
sx q[0];
rz(-1.3716475) q[0];
sx q[0];
rz(2.7535723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74036171) q[2];
sx q[2];
rz(-2.2380793) q[2];
sx q[2];
rz(2.0118535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.302205) q[1];
sx q[1];
rz(-2.6714462) q[1];
sx q[1];
rz(-0.66580982) q[1];
x q[2];
rz(0.42428942) q[3];
sx q[3];
rz(-2.9919713) q[3];
sx q[3];
rz(-0.77009088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55769354) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(1.0830797) q[2];
rz(-0.22100581) q[3];
sx q[3];
rz(-2.998304) q[3];
sx q[3];
rz(2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1204656) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(-0.13588151) q[0];
rz(-1.7120301) q[1];
sx q[1];
rz(-1.1212965) q[1];
sx q[1];
rz(2.8543465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3207239) q[0];
sx q[0];
rz(-2.831376) q[0];
sx q[0];
rz(-0.91281548) q[0];
rz(-1.9075228) q[2];
sx q[2];
rz(-0.25114608) q[2];
sx q[2];
rz(2.2876571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96480709) q[1];
sx q[1];
rz(-0.63399708) q[1];
sx q[1];
rz(0.37668677) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8494037) q[3];
sx q[3];
rz(-2.1351519) q[3];
sx q[3];
rz(-0.60114229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91934943) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(2.3075721) q[2];
rz(2.833994) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(-2.8789177) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.96497) q[0];
sx q[0];
rz(-1.8360092) q[0];
sx q[0];
rz(-1.0259884) q[0];
rz(0.38656825) q[1];
sx q[1];
rz(-1.760512) q[1];
sx q[1];
rz(-0.61813933) q[1];
rz(-1.7025399) q[2];
sx q[2];
rz(-1.9921039) q[2];
sx q[2];
rz(2.8767246) q[2];
rz(2.9959804) q[3];
sx q[3];
rz(-1.3619193) q[3];
sx q[3];
rz(-3.0633434) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
