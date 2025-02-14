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
rz(2.3020701) q[0];
sx q[0];
rz(-1.5538202) q[0];
sx q[0];
rz(2.177218) q[0];
rz(-2.7388465) q[1];
sx q[1];
rz(-1.4637113) q[1];
sx q[1];
rz(0.95199624) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6168686) q[0];
sx q[0];
rz(-1.476912) q[0];
sx q[0];
rz(2.062172) q[0];
rz(1.9929041) q[2];
sx q[2];
rz(-1.4808345) q[2];
sx q[2];
rz(0.46101704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9612572) q[1];
sx q[1];
rz(-1.9587659) q[1];
sx q[1];
rz(-1.9684102) q[1];
rz(0.49523109) q[3];
sx q[3];
rz(-0.81988813) q[3];
sx q[3];
rz(2.8748729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6160148) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(2.7363321) q[2];
rz(2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90527117) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(2.3765833) q[0];
rz(2.8969823) q[1];
sx q[1];
rz(-1.8965992) q[1];
sx q[1];
rz(0.6388706) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7981804) q[0];
sx q[0];
rz(-1.5454146) q[0];
sx q[0];
rz(2.4973386) q[0];
rz(2.6354796) q[2];
sx q[2];
rz(-0.68685907) q[2];
sx q[2];
rz(-3.0619301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7201405) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(-1.1790183) q[1];
rz(-1.8539393) q[3];
sx q[3];
rz(-2.3842709) q[3];
sx q[3];
rz(-2.9589096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8581533) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-2.8972674) q[2];
rz(-2.0325932) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(-2.9037156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7636488) q[0];
sx q[0];
rz(-1.0841333) q[0];
sx q[0];
rz(-1.2373244) q[0];
rz(1.411865) q[1];
sx q[1];
rz(-2.6928163) q[1];
sx q[1];
rz(1.3317187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004382172) q[0];
sx q[0];
rz(-1.6586275) q[0];
sx q[0];
rz(2.3669404) q[0];
rz(1.9218512) q[2];
sx q[2];
rz(-1.7760385) q[2];
sx q[2];
rz(2.308941) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0432555) q[1];
sx q[1];
rz(-1.6041293) q[1];
sx q[1];
rz(-0.055524932) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3490155) q[3];
sx q[3];
rz(-2.3537785) q[3];
sx q[3];
rz(0.28258309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1639013) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(-0.65995556) q[2];
rz(-1.4466977) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25809836) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(-2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(2.5925327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814007) q[0];
sx q[0];
rz(-2.9219919) q[0];
sx q[0];
rz(-0.85217969) q[0];
x q[1];
rz(1.8070776) q[2];
sx q[2];
rz(-2.5388814) q[2];
sx q[2];
rz(0.38379764) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.258865) q[1];
sx q[1];
rz(-2.8587864) q[1];
sx q[1];
rz(0.23587366) q[1];
x q[2];
rz(-1.9677343) q[3];
sx q[3];
rz(-1.8324495) q[3];
sx q[3];
rz(1.73571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.05222008) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(0.29903665) q[2];
rz(2.5236169) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(1.3581902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.95737902) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(-2.6755565) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(0.83494157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3059495) q[0];
sx q[0];
rz(-1.5993274) q[0];
sx q[0];
rz(3.0898407) q[0];
x q[1];
rz(-2.4136426) q[2];
sx q[2];
rz(-3.0228428) q[2];
sx q[2];
rz(0.70357067) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2920759) q[1];
sx q[1];
rz(-2.6585732) q[1];
sx q[1];
rz(1.6860875) q[1];
rz(-0.68170516) q[3];
sx q[3];
rz(-0.85452628) q[3];
sx q[3];
rz(-2.0264421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20778188) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(-0.031522838) q[2];
rz(0.75782123) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23444489) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(0.031540792) q[0];
rz(3.0517598) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(2.8057742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528362) q[0];
sx q[0];
rz(-1.2169219) q[0];
sx q[0];
rz(1.2214946) q[0];
x q[1];
rz(-1.01891) q[2];
sx q[2];
rz(-1.8621716) q[2];
sx q[2];
rz(-2.5138829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0586186) q[1];
sx q[1];
rz(-2.1095643) q[1];
sx q[1];
rz(-2.8584216) q[1];
rz(-1.0085513) q[3];
sx q[3];
rz(-2.6186071) q[3];
sx q[3];
rz(-2.6949835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.971656) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(-1.2706832) q[2];
rz(0.058567889) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(0.2659604) q[0];
rz(2.3174441) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(-1.5305653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8088246) q[0];
sx q[0];
rz(-0.36697373) q[0];
sx q[0];
rz(-0.76353188) q[0];
rz(-0.73101298) q[2];
sx q[2];
rz(-1.5314529) q[2];
sx q[2];
rz(0.61737663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8083473) q[1];
sx q[1];
rz(-2.5032804) q[1];
sx q[1];
rz(0.45175938) q[1];
x q[2];
rz(2.7830151) q[3];
sx q[3];
rz(-2.3365575) q[3];
sx q[3];
rz(0.81881675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.019235762) q[2];
sx q[2];
rz(-2.0255721) q[2];
sx q[2];
rz(2.1717211) q[2];
rz(-2.9979749) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(-0.6374878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.403991) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(-0.10424374) q[0];
rz(-0.741611) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(-0.43824497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40566984) q[0];
sx q[0];
rz(-0.68016648) q[0];
sx q[0];
rz(-0.16323547) q[0];
x q[1];
rz(-2.1684219) q[2];
sx q[2];
rz(-1.3750374) q[2];
sx q[2];
rz(-2.6225066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91854535) q[1];
sx q[1];
rz(-0.96541903) q[1];
sx q[1];
rz(0.041408509) q[1];
x q[2];
rz(-1.7167818) q[3];
sx q[3];
rz(-2.101482) q[3];
sx q[3];
rz(2.3455181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7354342) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(0.3212277) q[2];
rz(0.99212232) q[3];
sx q[3];
rz(-1.0935874) q[3];
sx q[3];
rz(1.4120302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883437) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(-0.14800063) q[0];
rz(-2.0026228) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(0.99050561) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73511998) q[0];
sx q[0];
rz(-0.48790259) q[0];
sx q[0];
rz(0.16012971) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6338088) q[2];
sx q[2];
rz(-2.3427026) q[2];
sx q[2];
rz(0.2474167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3504178) q[1];
sx q[1];
rz(-1.5656078) q[1];
sx q[1];
rz(-1.6020011) q[1];
rz(-pi) q[2];
rz(2.7630291) q[3];
sx q[3];
rz(-1.1051911) q[3];
sx q[3];
rz(-1.9052802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(-2.5992744) q[2];
rz(-1.6769241) q[3];
sx q[3];
rz(-1.9138391) q[3];
sx q[3];
rz(0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2391613) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-2.8620128) q[0];
rz(-2.3565049) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(-0.40207544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063906625) q[0];
sx q[0];
rz(-0.55372596) q[0];
sx q[0];
rz(2.5095449) q[0];
x q[1];
rz(0.80569141) q[2];
sx q[2];
rz(-2.1955954) q[2];
sx q[2];
rz(0.8936409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6732503) q[1];
sx q[1];
rz(-1.3788917) q[1];
sx q[1];
rz(0.90190355) q[1];
rz(0.28487408) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(-1.3314825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98099199) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(-0.92657363) q[2];
rz(2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.8615104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67113039) q[0];
sx q[0];
rz(-2.0736546) q[0];
sx q[0];
rz(0.14508844) q[0];
rz(0.078770272) q[1];
sx q[1];
rz(-1.7369743) q[1];
sx q[1];
rz(-1.8440934) q[1];
rz(-3.1263951) q[2];
sx q[2];
rz(-0.79833818) q[2];
sx q[2];
rz(-1.0865059) q[2];
rz(1.3435626) q[3];
sx q[3];
rz(-0.5323635) q[3];
sx q[3];
rz(-1.3743286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
