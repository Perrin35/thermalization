OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1998347) q[0];
sx q[0];
rz(-0.69154843) q[0];
sx q[0];
rz(-0.77342311) q[0];
rz(3.050488) q[1];
sx q[1];
rz(2.3478822) q[1];
sx q[1];
rz(4.4749727) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0036732) q[0];
sx q[0];
rz(-0.9124476) q[0];
sx q[0];
rz(-1.5821032) q[0];
rz(-pi) q[1];
rz(-2.0699507) q[2];
sx q[2];
rz(-1.3678275) q[2];
sx q[2];
rz(0.90107337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7385834) q[1];
sx q[1];
rz(-2.3883551) q[1];
sx q[1];
rz(2.1639216) q[1];
rz(-pi) q[2];
rz(3.1245272) q[3];
sx q[3];
rz(-1.6380042) q[3];
sx q[3];
rz(-0.9931356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8495463) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(0.41479659) q[2];
rz(-1.8997806) q[3];
sx q[3];
rz(-1.1153699) q[3];
sx q[3];
rz(-0.19150664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5135797) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(0.43981788) q[0];
rz(-3.0385333) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(-1.1691079) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1531853) q[0];
sx q[0];
rz(-2.4988032) q[0];
sx q[0];
rz(-1.8039186) q[0];
x q[1];
rz(-2.1610519) q[2];
sx q[2];
rz(-1.613918) q[2];
sx q[2];
rz(1.0999964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1406469) q[1];
sx q[1];
rz(-0.92921153) q[1];
sx q[1];
rz(-2.916275) q[1];
rz(-pi) q[2];
rz(0.86128791) q[3];
sx q[3];
rz(-1.420974) q[3];
sx q[3];
rz(3.0784208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.64572) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(0.24620852) q[2];
rz(-1.5919033) q[3];
sx q[3];
rz(-0.61384765) q[3];
sx q[3];
rz(-1.3932108) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.15568) q[0];
sx q[0];
rz(-2.0014626) q[0];
sx q[0];
rz(2.929856) q[0];
rz(3.1302997) q[1];
sx q[1];
rz(-2.3432422) q[1];
sx q[1];
rz(1.4439772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31682184) q[0];
sx q[0];
rz(-1.5708692) q[0];
sx q[0];
rz(-3.1064338) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3668187) q[2];
sx q[2];
rz(-2.8958671) q[2];
sx q[2];
rz(-1.7468921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8577598) q[1];
sx q[1];
rz(-0.74537828) q[1];
sx q[1];
rz(0.56826313) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6908247) q[3];
sx q[3];
rz(-1.4891461) q[3];
sx q[3];
rz(-1.1686366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6649365) q[2];
sx q[2];
rz(-1.0756476) q[2];
sx q[2];
rz(-2.5962489) q[2];
rz(0.90965811) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(1.9285412) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.720541) q[0];
sx q[0];
rz(-2.4626829) q[0];
sx q[0];
rz(-2.3140267) q[0];
rz(-1.958581) q[1];
sx q[1];
rz(-1.8667826) q[1];
sx q[1];
rz(-1.8756728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1249378) q[0];
sx q[0];
rz(-1.6724574) q[0];
sx q[0];
rz(0.09135017) q[0];
x q[1];
rz(-1.363344) q[2];
sx q[2];
rz(-1.5816763) q[2];
sx q[2];
rz(1.915773) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79252473) q[1];
sx q[1];
rz(-1.6025839) q[1];
sx q[1];
rz(2.1963901) q[1];
rz(-pi) q[2];
rz(-1.7608579) q[3];
sx q[3];
rz(-1.7281088) q[3];
sx q[3];
rz(-0.33657956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2563235) q[2];
sx q[2];
rz(-1.6689391) q[2];
sx q[2];
rz(3.0598158) q[2];
rz(-2.8335588) q[3];
sx q[3];
rz(-0.48397288) q[3];
sx q[3];
rz(1.6035732) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643739) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(0.070505738) q[0];
rz(0.33440691) q[1];
sx q[1];
rz(-0.53135482) q[1];
sx q[1];
rz(-2.6883584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.189126) q[0];
sx q[0];
rz(-1.0116018) q[0];
sx q[0];
rz(0.97049148) q[0];
rz(-2.9137827) q[2];
sx q[2];
rz(-1.1880298) q[2];
sx q[2];
rz(-1.5441976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91562781) q[1];
sx q[1];
rz(-2.8702822) q[1];
sx q[1];
rz(2.7167999) q[1];
x q[2];
rz(-2.2599561) q[3];
sx q[3];
rz(-2.305718) q[3];
sx q[3];
rz(0.78208941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8411023) q[2];
sx q[2];
rz(-1.5056242) q[2];
sx q[2];
rz(0.21031586) q[2];
rz(2.7503843) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(-0.44262639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5000358) q[0];
sx q[0];
rz(-1.1316725) q[0];
sx q[0];
rz(0.17182194) q[0];
rz(-1.7956644) q[1];
sx q[1];
rz(-2.185952) q[1];
sx q[1];
rz(1.1489493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0003819) q[0];
sx q[0];
rz(-0.11185574) q[0];
sx q[0];
rz(3.0684708) q[0];
x q[1];
rz(0.8900155) q[2];
sx q[2];
rz(-0.35951722) q[2];
sx q[2];
rz(0.36717626) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90075452) q[1];
sx q[1];
rz(-0.60610702) q[1];
sx q[1];
rz(-1.89487) q[1];
rz(-pi) q[2];
rz(-0.20307417) q[3];
sx q[3];
rz(-1.0725029) q[3];
sx q[3];
rz(-0.6607252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12080869) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(0.42050335) q[2];
rz(-1.486982) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(1.9871064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7789975) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(-2.763789) q[0];
rz(1.3899577) q[1];
sx q[1];
rz(-1.708834) q[1];
sx q[1];
rz(-0.43209824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2588151) q[0];
sx q[0];
rz(-1.3963789) q[0];
sx q[0];
rz(2.9091878) q[0];
x q[1];
rz(-0.092707002) q[2];
sx q[2];
rz(-1.6039298) q[2];
sx q[2];
rz(-0.2011782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3258265) q[1];
sx q[1];
rz(-1.5349421) q[1];
sx q[1];
rz(0.46200606) q[1];
rz(-pi) q[2];
rz(2.2446052) q[3];
sx q[3];
rz(-1.8963061) q[3];
sx q[3];
rz(1.432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3720588) q[2];
sx q[2];
rz(-2.5199315) q[2];
sx q[2];
rz(0.2335693) q[2];
rz(-1.6893049) q[3];
sx q[3];
rz(-1.4315616) q[3];
sx q[3];
rz(1.5610032) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643144) q[0];
sx q[0];
rz(-1.0294585) q[0];
sx q[0];
rz(2.8218063) q[0];
rz(-0.86209595) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(1.7822942) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0727973) q[0];
sx q[0];
rz(-0.85155838) q[0];
sx q[0];
rz(-0.79049514) q[0];
rz(0.14131693) q[2];
sx q[2];
rz(-1.3719146) q[2];
sx q[2];
rz(0.74706739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2923802) q[1];
sx q[1];
rz(-1.74382) q[1];
sx q[1];
rz(-2.4166957) q[1];
x q[2];
rz(0.93710758) q[3];
sx q[3];
rz(-1.345063) q[3];
sx q[3];
rz(-2.2645775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84997815) q[2];
sx q[2];
rz(-0.41899592) q[2];
sx q[2];
rz(-2.6739547) q[2];
rz(-0.48401287) q[3];
sx q[3];
rz(-1.7734807) q[3];
sx q[3];
rz(1.2776933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17396946) q[0];
sx q[0];
rz(-1.6521709) q[0];
sx q[0];
rz(-2.0066579) q[0];
rz(-3.0813772) q[1];
sx q[1];
rz(-2.4048012) q[1];
sx q[1];
rz(-1.5312451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519015) q[0];
sx q[0];
rz(-0.78329059) q[0];
sx q[0];
rz(2.2499311) q[0];
rz(-pi) q[1];
rz(2.5217053) q[2];
sx q[2];
rz(-2.3522178) q[2];
sx q[2];
rz(0.7482341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.798315) q[1];
sx q[1];
rz(-0.28099842) q[1];
sx q[1];
rz(2.7839989) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6100735) q[3];
sx q[3];
rz(-0.78900064) q[3];
sx q[3];
rz(1.5189182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81426364) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(0.72861707) q[2];
rz(0.21555756) q[3];
sx q[3];
rz(-1.39648) q[3];
sx q[3];
rz(-1.0936776) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88290596) q[0];
sx q[0];
rz(-0.031143324) q[0];
sx q[0];
rz(-2.8236142) q[0];
rz(-1.7440354) q[1];
sx q[1];
rz(-1.3293068) q[1];
sx q[1];
rz(-0.2831645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4508789) q[0];
sx q[0];
rz(-0.15199272) q[0];
sx q[0];
rz(-0.99894036) q[0];
rz(-0.0816143) q[2];
sx q[2];
rz(-0.75576648) q[2];
sx q[2];
rz(2.4163742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8258543) q[1];
sx q[1];
rz(-0.87402499) q[1];
sx q[1];
rz(2.524074) q[1];
rz(-pi) q[2];
rz(3.0491676) q[3];
sx q[3];
rz(-1.7778977) q[3];
sx q[3];
rz(0.3849808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2240923) q[2];
sx q[2];
rz(-1.6250236) q[2];
sx q[2];
rz(0.089281233) q[2];
rz(1.3034405) q[3];
sx q[3];
rz(-2.3035514) q[3];
sx q[3];
rz(2.1681521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3598809) q[0];
sx q[0];
rz(-0.64890535) q[0];
sx q[0];
rz(-2.1606408) q[0];
rz(-2.5557062) q[1];
sx q[1];
rz(-1.2955019) q[1];
sx q[1];
rz(-1.6265709) q[1];
rz(0.21296756) q[2];
sx q[2];
rz(-1.7473379) q[2];
sx q[2];
rz(2.9416549) q[2];
rz(1.1896776) q[3];
sx q[3];
rz(-1.3314691) q[3];
sx q[3];
rz(-0.21503147) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
