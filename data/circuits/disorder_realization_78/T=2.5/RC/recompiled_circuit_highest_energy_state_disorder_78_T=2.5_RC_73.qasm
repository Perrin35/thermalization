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
rz(-3.0566293) q[0];
sx q[0];
rz(-0.30240107) q[0];
sx q[0];
rz(-3.1095355) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(3.6454522) q[1];
sx q[1];
rz(7.03581) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0467211) q[0];
sx q[0];
rz(-1.7716452) q[0];
sx q[0];
rz(-0.3263536) q[0];
x q[1];
rz(-1.1967802) q[2];
sx q[2];
rz(-1.2311282) q[2];
sx q[2];
rz(2.938478) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51576534) q[1];
sx q[1];
rz(-1.2316362) q[1];
sx q[1];
rz(1.8372507) q[1];
rz(-pi) q[2];
rz(2.9273622) q[3];
sx q[3];
rz(-1.9211244) q[3];
sx q[3];
rz(-0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(0.44102937) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(-0.41020694) q[0];
rz(-0.38093105) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(1.7832696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5986431) q[0];
sx q[0];
rz(-0.99587593) q[0];
sx q[0];
rz(-0.47358124) q[0];
rz(-pi) q[1];
rz(2.3097281) q[2];
sx q[2];
rz(-0.68507776) q[2];
sx q[2];
rz(-2.5545504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4315003) q[1];
sx q[1];
rz(-2.9572801) q[1];
sx q[1];
rz(1.1460365) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3639327) q[3];
sx q[3];
rz(-2.4914722) q[3];
sx q[3];
rz(-0.93322414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(0.50296339) q[2];
rz(-1.2991692) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-0.51173156) q[0];
sx q[0];
rz(0.9758392) q[0];
rz(2.2052235) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(-0.13793129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9741623) q[0];
sx q[0];
rz(-0.61744055) q[0];
sx q[0];
rz(-2.2327044) q[0];
rz(1.6866272) q[2];
sx q[2];
rz(-1.9660913) q[2];
sx q[2];
rz(2.9651053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5441078) q[1];
sx q[1];
rz(-2.4307837) q[1];
sx q[1];
rz(1.3686468) q[1];
rz(1.2983599) q[3];
sx q[3];
rz(-2.0217293) q[3];
sx q[3];
rz(-0.48905269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5522573) q[2];
sx q[2];
rz(-1.5312803) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(0.79834437) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524356) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(2.074746) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(-2.7222395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28946149) q[0];
sx q[0];
rz(-1.7999819) q[0];
sx q[0];
rz(0.16332345) q[0];
rz(-pi) q[1];
rz(-1.4801976) q[2];
sx q[2];
rz(-1.0444469) q[2];
sx q[2];
rz(-1.3535997) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90430656) q[1];
sx q[1];
rz(-1.6349873) q[1];
sx q[1];
rz(-1.2130846) q[1];
rz(-2.0556424) q[3];
sx q[3];
rz(-1.2618229) q[3];
sx q[3];
rz(-3.0687208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7833917) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(-1.5436714) q[2];
rz(-1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.876038) q[0];
sx q[0];
rz(-0.72172481) q[0];
sx q[0];
rz(-0.75620404) q[0];
rz(1.8887695) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(-1.7255712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21391103) q[0];
sx q[0];
rz(-2.1577303) q[0];
sx q[0];
rz(-1.3089433) q[0];
x q[1];
rz(0.33966392) q[2];
sx q[2];
rz(-0.36298266) q[2];
sx q[2];
rz(-1.5756922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6452613) q[1];
sx q[1];
rz(-1.0320391) q[1];
sx q[1];
rz(-1.5625728) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6719499) q[3];
sx q[3];
rz(-2.7386463) q[3];
sx q[3];
rz(-2.6279891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8004134) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(2.9608012) q[2];
rz(-1.6527269) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(-1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(-0.80023009) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-3.0175993) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127045) q[0];
sx q[0];
rz(-1.4409954) q[0];
sx q[0];
rz(-0.028667269) q[0];
rz(-3.0818865) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(1.6725181) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2468894) q[1];
sx q[1];
rz(-1.2661625) q[1];
sx q[1];
rz(1.4851634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97449755) q[3];
sx q[3];
rz(-1.8829405) q[3];
sx q[3];
rz(-2.6016935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2064712) q[2];
sx q[2];
rz(-1.7266885) q[2];
sx q[2];
rz(0.07621152) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(0.61856234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.16977075) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(2.3845657) q[0];
rz(0.79611671) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(-0.98181358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1025006) q[0];
sx q[0];
rz(-1.8538686) q[0];
sx q[0];
rz(-2.3209178) q[0];
rz(-2.1481291) q[2];
sx q[2];
rz(-3.0138123) q[2];
sx q[2];
rz(1.8005231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0701778) q[1];
sx q[1];
rz(-1.1648253) q[1];
sx q[1];
rz(2.0858048) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8514093) q[3];
sx q[3];
rz(-1.4320489) q[3];
sx q[3];
rz(-0.25158238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7299399) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(-2.5977503) q[2];
rz(0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.2127317) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(-2.5182356) q[0];
rz(-0.32132545) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(2.8772433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0118222) q[0];
sx q[0];
rz(-1.3589824) q[0];
sx q[0];
rz(2.2925633) q[0];
rz(-pi) q[1];
rz(0.041548403) q[2];
sx q[2];
rz(-1.9061521) q[2];
sx q[2];
rz(-1.2183508) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.2920989) q[1];
sx q[1];
rz(2.3257117) q[1];
rz(-3.0086413) q[3];
sx q[3];
rz(-1.2556517) q[3];
sx q[3];
rz(-0.12115762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0048206) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(0.36561203) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356165) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-2.0756508) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(-0.33445439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850824) q[0];
sx q[0];
rz(-1.2691536) q[0];
sx q[0];
rz(-1.3972052) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43844079) q[2];
sx q[2];
rz(-1.6980357) q[2];
sx q[2];
rz(-1.8220937) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0171256) q[1];
sx q[1];
rz(-1.4216058) q[1];
sx q[1];
rz(-2.6890432) q[1];
rz(0.79631898) q[3];
sx q[3];
rz(-1.2346754) q[3];
sx q[3];
rz(1.6567865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(3.0706792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(2.8712414) q[0];
rz(-1.5125754) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975166) q[0];
sx q[0];
rz(-2.368921) q[0];
sx q[0];
rz(-2.0324043) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43105189) q[2];
sx q[2];
rz(-2.1630175) q[2];
sx q[2];
rz(2.2811802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8737265) q[1];
sx q[1];
rz(-1.6609816) q[1];
sx q[1];
rz(-1.6117881) q[1];
rz(-0.38624951) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(-1.7084427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35499972) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(-2.2641838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56562051) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(1.3427973) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-3.1411896) q[2];
sx q[2];
rz(-0.12231356) q[2];
sx q[2];
rz(-0.077153645) q[2];
rz(-2.9407637) q[3];
sx q[3];
rz(-2.8786537) q[3];
sx q[3];
rz(-2.636551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
