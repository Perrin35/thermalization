OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(-2.3242216) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(-1.2004381) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49657208) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(2.7215331) q[0];
rz(-pi) q[1];
rz(1.0317694) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(-0.10345085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8564954) q[1];
sx q[1];
rz(-0.76438099) q[1];
sx q[1];
rz(-2.7544423) q[1];
x q[2];
rz(-1.2458385) q[3];
sx q[3];
rz(-2.339139) q[3];
sx q[3];
rz(1.6368395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.583064) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(3.0875207) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(-0.53584677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(2.9931086) q[0];
x q[1];
rz(1.3185805) q[2];
sx q[2];
rz(-2.2644342) q[2];
sx q[2];
rz(-2.0351621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9745001) q[1];
sx q[1];
rz(-1.316861) q[1];
sx q[1];
rz(-0.42029917) q[1];
rz(1.4278533) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-0.95834857) q[2];
rz(-3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31608554) q[0];
sx q[0];
rz(-1.320991) q[0];
sx q[0];
rz(-0.020629701) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0518968) q[2];
sx q[2];
rz(-1.2604453) q[2];
sx q[2];
rz(0.87583625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3734696) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(2.5440689) q[1];
rz(0.13173007) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4767544) q[0];
sx q[0];
rz(-1.1676482) q[0];
sx q[0];
rz(-0.82157764) q[0];
rz(2.9085607) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(-0.88569966) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.199898) q[1];
sx q[1];
rz(-2.3823793) q[1];
sx q[1];
rz(-0.56337507) q[1];
x q[2];
rz(-0.36066182) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(0.053645596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-1.0446576) q[2];
rz(-0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(-0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9263822) q[0];
sx q[0];
rz(-0.97391093) q[0];
sx q[0];
rz(-2.9125288) q[0];
rz(0.28508614) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(-1.8269055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16380285) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(2.8603641) q[1];
x q[2];
rz(-0.78955663) q[3];
sx q[3];
rz(-1.4649179) q[3];
sx q[3];
rz(-1.4279799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832702) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(-1.8002611) q[0];
rz(-pi) q[1];
rz(2.5541359) q[2];
sx q[2];
rz(-0.60411462) q[2];
sx q[2];
rz(0.47746745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7212972) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(0.23423127) q[1];
rz(-pi) q[2];
rz(-0.22535725) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(-0.55707896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6713312) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11628843) q[0];
sx q[0];
rz(-1.3447273) q[0];
sx q[0];
rz(-0.53209214) q[0];
x q[1];
rz(0.62692554) q[2];
sx q[2];
rz(-1.727384) q[2];
sx q[2];
rz(2.2121034) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6262498) q[1];
sx q[1];
rz(-1.5222349) q[1];
sx q[1];
rz(1.6471144) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-0.75546414) q[3];
sx q[3];
rz(-1.2681703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1091052) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(-2.0672654) q[0];
rz(1.9365963) q[2];
sx q[2];
rz(-1.6779643) q[2];
sx q[2];
rz(-2.0049713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3216746) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(0.8701156) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91248625) q[3];
sx q[3];
rz(-2.6819326) q[3];
sx q[3];
rz(-2.130079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(-0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-2.2414482) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(0.94690698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3627975) q[0];
sx q[0];
rz(-2.6914094) q[0];
sx q[0];
rz(2.3953526) q[0];
x q[1];
rz(-0.38893716) q[2];
sx q[2];
rz(-0.10483327) q[2];
sx q[2];
rz(-0.36738415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.088989181) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(-2.8903972) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9506748) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-0.65504909) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(-2.8331579) q[0];
x q[1];
rz(-1.3078493) q[2];
sx q[2];
rz(-0.73336468) q[2];
sx q[2];
rz(-2.3761689) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6380784) q[1];
sx q[1];
rz(-0.44736171) q[1];
sx q[1];
rz(-0.7034941) q[1];
rz(-pi) q[2];
x q[2];
rz(2.321407) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(1.5370777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(2.9121493) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-1.3953801) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(2.2453528) q[3];
sx q[3];
rz(-1.8825718) q[3];
sx q[3];
rz(2.752302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];