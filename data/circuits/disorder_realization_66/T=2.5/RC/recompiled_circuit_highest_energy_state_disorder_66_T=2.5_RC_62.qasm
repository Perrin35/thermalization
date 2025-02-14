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
rz(-1.00296) q[0];
sx q[0];
rz(-2.668219) q[0];
sx q[0];
rz(1.0869429) q[0];
rz(-0.76397693) q[1];
sx q[1];
rz(-0.44372258) q[1];
sx q[1];
rz(-0.8313764) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87693857) q[0];
sx q[0];
rz(-2.0994517) q[0];
sx q[0];
rz(-2.1318046) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23687266) q[2];
sx q[2];
rz(-1.1129967) q[2];
sx q[2];
rz(-1.906369) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8352141) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(1.4500965) q[1];
x q[2];
rz(3.0696297) q[3];
sx q[3];
rz(-1.5441893) q[3];
sx q[3];
rz(0.047544971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5272556) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(0.41198507) q[2];
rz(0.24474239) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(-2.8542724) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7016542) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(2.7213851) q[0];
rz(-0.48928753) q[1];
sx q[1];
rz(-0.91541618) q[1];
sx q[1];
rz(1.9583826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1294043) q[0];
sx q[0];
rz(-0.79372915) q[0];
sx q[0];
rz(0.94214006) q[0];
x q[1];
rz(2.8274306) q[2];
sx q[2];
rz(-1.2622241) q[2];
sx q[2];
rz(-1.6752288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.015731363) q[1];
sx q[1];
rz(-2.1441048) q[1];
sx q[1];
rz(2.0290976) q[1];
x q[2];
rz(1.7443827) q[3];
sx q[3];
rz(-1.6866444) q[3];
sx q[3];
rz(2.4811452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6861787) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(-0.14032826) q[2];
rz(2.8651107) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(-2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8088733) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(-1.3746388) q[0];
rz(-2.2287492) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(2.1106145) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46488097) q[0];
sx q[0];
rz(-1.2455938) q[0];
sx q[0];
rz(-2.9487378) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1616251) q[2];
sx q[2];
rz(-2.4942932) q[2];
sx q[2];
rz(1.5673296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7056641) q[1];
sx q[1];
rz(-1.9829653) q[1];
sx q[1];
rz(-0.010910587) q[1];
rz(2.616373) q[3];
sx q[3];
rz(-0.37074836) q[3];
sx q[3];
rz(1.4300089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-2.5035109) q[2];
rz(-3.0408527) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(0.49873275) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8345555) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(1.3822973) q[0];
rz(2.8221829) q[1];
sx q[1];
rz(-1.6326135) q[1];
sx q[1];
rz(-2.3841948) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7865407) q[0];
sx q[0];
rz(-3.0664223) q[0];
sx q[0];
rz(3.0274903) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73586936) q[2];
sx q[2];
rz(-0.19093787) q[2];
sx q[2];
rz(1.5895933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0075052) q[1];
sx q[1];
rz(-2.0148835) q[1];
sx q[1];
rz(0.92625463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0946461) q[3];
sx q[3];
rz(-2.110384) q[3];
sx q[3];
rz(1.7612518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77183977) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(-2.3337951) q[2];
rz(-1.7317023) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(-0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(2.156303) q[0];
rz(1.4957042) q[1];
sx q[1];
rz(-2.0031877) q[1];
sx q[1];
rz(-2.2122502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083223) q[0];
sx q[0];
rz(-2.4886971) q[0];
sx q[0];
rz(0.44954957) q[0];
x q[1];
rz(-2.5856951) q[2];
sx q[2];
rz(-2.4553442) q[2];
sx q[2];
rz(-2.2874119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29581901) q[1];
sx q[1];
rz(-2.5717178) q[1];
sx q[1];
rz(-1.946158) q[1];
rz(-1.0096022) q[3];
sx q[3];
rz(-1.582904) q[3];
sx q[3];
rz(-0.59134134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0890961) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(2.055114) q[2];
rz(0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(1.121678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1679967) q[0];
sx q[0];
rz(-0.43742988) q[0];
sx q[0];
rz(-0.22015372) q[0];
rz(1.5036229) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(-0.55353177) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93388825) q[0];
sx q[0];
rz(-1.7676864) q[0];
sx q[0];
rz(-1.4536727) q[0];
rz(-pi) q[1];
rz(0.4398063) q[2];
sx q[2];
rz(-2.7436119) q[2];
sx q[2];
rz(2.8619638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11343918) q[1];
sx q[1];
rz(-2.4467797) q[1];
sx q[1];
rz(-3.0656612) q[1];
rz(-pi) q[2];
rz(-2.8643478) q[3];
sx q[3];
rz(-2.436815) q[3];
sx q[3];
rz(0.98990209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0452051) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(0.20827797) q[2];
rz(1.1194718) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(-2.313405) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49523062) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(1.0736504) q[0];
rz(-0.13058361) q[1];
sx q[1];
rz(-1.5740296) q[1];
sx q[1];
rz(-2.1317587) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829941) q[0];
sx q[0];
rz(-1.3855591) q[0];
sx q[0];
rz(0.40249975) q[0];
x q[1];
rz(-2.9931941) q[2];
sx q[2];
rz(-1.6434533) q[2];
sx q[2];
rz(2.1138482) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9562137) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(0.57422178) q[1];
rz(-pi) q[2];
rz(-2.4773682) q[3];
sx q[3];
rz(-1.281732) q[3];
sx q[3];
rz(-2.5632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4449571) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(-0.32336393) q[2];
rz(-2.4671386) q[3];
sx q[3];
rz(-0.89265299) q[3];
sx q[3];
rz(-3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12482878) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-2.0523409) q[0];
rz(2.543154) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(-2.3425897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4161637) q[0];
sx q[0];
rz(-1.874864) q[0];
sx q[0];
rz(1.2771525) q[0];
rz(1.4510462) q[2];
sx q[2];
rz(-0.66968067) q[2];
sx q[2];
rz(1.0543752) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1027229) q[1];
sx q[1];
rz(-0.33867044) q[1];
sx q[1];
rz(-2.6790957) q[1];
x q[2];
rz(0.022074583) q[3];
sx q[3];
rz(-0.24188731) q[3];
sx q[3];
rz(2.8245408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4716855) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(-1.8890107) q[2];
rz(-2.0910828) q[3];
sx q[3];
rz(-0.59058467) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893386) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(-2.6259212) q[0];
rz(-1.6853261) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(-2.8175443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0312248) q[0];
sx q[0];
rz(-0.47981167) q[0];
sx q[0];
rz(-2.5672002) q[0];
x q[1];
rz(0.44948795) q[2];
sx q[2];
rz(-0.96908334) q[2];
sx q[2];
rz(0.692761) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1157323) q[1];
sx q[1];
rz(-2.6792791) q[1];
sx q[1];
rz(1.8883392) q[1];
rz(-2.1284655) q[3];
sx q[3];
rz(-1.6383759) q[3];
sx q[3];
rz(2.9014719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9825762) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(2.4328361) q[2];
rz(-0.7835663) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.5172575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.5358955) q[0];
sx q[0];
rz(-0.63740969) q[0];
sx q[0];
rz(-0.15644431) q[0];
rz(0.63977301) q[1];
sx q[1];
rz(-2.4867058) q[1];
sx q[1];
rz(0.99008647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759686) q[0];
sx q[0];
rz(-0.98593119) q[0];
sx q[0];
rz(0.82230277) q[0];
rz(-pi) q[1];
rz(-1.9205356) q[2];
sx q[2];
rz(-1.4863401) q[2];
sx q[2];
rz(-1.8366369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2459979) q[1];
sx q[1];
rz(-0.57064547) q[1];
sx q[1];
rz(2.8938608) q[1];
rz(-0.20179468) q[3];
sx q[3];
rz(-1.4425689) q[3];
sx q[3];
rz(1.1819672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.4116776) q[2];
sx q[2];
rz(-0.61335316) q[2];
rz(-2.8776045) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-0.67155513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(1.175066) q[2];
sx q[2];
rz(-1.9446951) q[2];
sx q[2];
rz(-2.4315628) q[2];
rz(2.522884) q[3];
sx q[3];
rz(-1.6029463) q[3];
sx q[3];
rz(-1.1986986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
