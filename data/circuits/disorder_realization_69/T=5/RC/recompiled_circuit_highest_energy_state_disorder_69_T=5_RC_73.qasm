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
rz(0.34613553) q[0];
sx q[0];
rz(-1.5032285) q[0];
sx q[0];
rz(0.7315973) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(-0.43391689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91125488) q[0];
sx q[0];
rz(-2.3106835) q[0];
sx q[0];
rz(-0.047538443) q[0];
rz(-pi) q[1];
rz(-2.5662759) q[2];
sx q[2];
rz(-1.0566718) q[2];
sx q[2];
rz(0.34808394) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9508236) q[1];
sx q[1];
rz(-2.1956823) q[1];
sx q[1];
rz(-2.6527219) q[1];
x q[2];
rz(-1.3796666) q[3];
sx q[3];
rz(-1.5401296) q[3];
sx q[3];
rz(-0.14740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71243858) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(2.1377371) q[2];
rz(-2.6058274) q[3];
sx q[3];
rz(-0.86447132) q[3];
sx q[3];
rz(0.0011477688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4681604) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(2.4144507) q[0];
rz(-0.06761059) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(-1.1236069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0824995) q[0];
sx q[0];
rz(-1.8525665) q[0];
sx q[0];
rz(2.6450794) q[0];
rz(-pi) q[1];
rz(-0.065326377) q[2];
sx q[2];
rz(-1.4484754) q[2];
sx q[2];
rz(2.9768012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.20642463) q[1];
sx q[1];
rz(-1.4713396) q[1];
sx q[1];
rz(1.8518492) q[1];
rz(-0.43895841) q[3];
sx q[3];
rz(-1.135313) q[3];
sx q[3];
rz(1.8451231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2850538) q[2];
sx q[2];
rz(-1.5679789) q[2];
sx q[2];
rz(2.6598568) q[2];
rz(-2.5887865) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8168617) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(3.0480296) q[0];
rz(-1.7612673) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(-2.5043452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940612) q[0];
sx q[0];
rz(-1.1406745) q[0];
sx q[0];
rz(-2.568666) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3445417) q[2];
sx q[2];
rz(-2.166894) q[2];
sx q[2];
rz(1.4486194) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6031076) q[1];
sx q[1];
rz(-2.7432495) q[1];
sx q[1];
rz(1.6377691) q[1];
x q[2];
rz(-1.5277385) q[3];
sx q[3];
rz(-1.2803439) q[3];
sx q[3];
rz(-1.8990979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20643413) q[2];
sx q[2];
rz(-2.5863402) q[2];
sx q[2];
rz(2.4578102) q[2];
rz(2.911496) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(2.7195215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049659599) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(-1.7012713) q[0];
rz(2.6662042) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-0.58116523) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6027045) q[0];
sx q[0];
rz(-1.4853444) q[0];
sx q[0];
rz(-2.2936725) q[0];
rz(-2.0330795) q[2];
sx q[2];
rz(-1.7812742) q[2];
sx q[2];
rz(-2.0049948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46266471) q[1];
sx q[1];
rz(-1.8724672) q[1];
sx q[1];
rz(0.48120166) q[1];
x q[2];
rz(2.851296) q[3];
sx q[3];
rz(-0.50945849) q[3];
sx q[3];
rz(-0.51028937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0873969) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(0.52097121) q[2];
rz(1.6446796) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(-1.6269256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5533376) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(2.0113373) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(-2.030453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53752995) q[0];
sx q[0];
rz(-1.1141234) q[0];
sx q[0];
rz(0.18484331) q[0];
rz(-pi) q[1];
rz(3.0658998) q[2];
sx q[2];
rz(-2.6106129) q[2];
sx q[2];
rz(0.39226433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9359253) q[1];
sx q[1];
rz(-1.7280518) q[1];
sx q[1];
rz(-2.2155511) q[1];
rz(-pi) q[2];
rz(-0.91577282) q[3];
sx q[3];
rz(-2.0213599) q[3];
sx q[3];
rz(2.5289867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7378716) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(-2.1007288) q[2];
rz(-0.8199842) q[3];
sx q[3];
rz(-1.2330202) q[3];
sx q[3];
rz(-2.5177054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(0.65648055) q[0];
rz(1.3735324) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(1.1318644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362909) q[0];
sx q[0];
rz(-1.3735008) q[0];
sx q[0];
rz(2.4407766) q[0];
rz(-2.4651338) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(1.187834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1161728) q[1];
sx q[1];
rz(-2.4301404) q[1];
sx q[1];
rz(1.6132212) q[1];
rz(-1.5664212) q[3];
sx q[3];
rz(-1.4203912) q[3];
sx q[3];
rz(0.37585092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6953096) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(-0.31965762) q[2];
rz(2.9122635) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(0.91013175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0685773) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(1.9101494) q[0];
rz(-3.0629509) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(-2.9873649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20040266) q[0];
sx q[0];
rz(-1.909167) q[0];
sx q[0];
rz(1.75941) q[0];
x q[1];
rz(1.9811822) q[2];
sx q[2];
rz(-1.850381) q[2];
sx q[2];
rz(-1.9094163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3992622) q[1];
sx q[1];
rz(-0.62472099) q[1];
sx q[1];
rz(-0.66466753) q[1];
x q[2];
rz(-0.46040543) q[3];
sx q[3];
rz(-0.45259991) q[3];
sx q[3];
rz(2.3456338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8057033) q[2];
sx q[2];
rz(-1.5259537) q[2];
sx q[2];
rz(-1.6501144) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(1.570805) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989646) q[0];
sx q[0];
rz(-2.0476116) q[0];
sx q[0];
rz(1.6866823) q[0];
rz(-0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(1.3386493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1478302) q[0];
sx q[0];
rz(-2.2394132) q[0];
sx q[0];
rz(2.6385959) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0089802) q[2];
sx q[2];
rz(-0.57719165) q[2];
sx q[2];
rz(-1.4711589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6852831) q[1];
sx q[1];
rz(-1.7830308) q[1];
sx q[1];
rz(2.2077399) q[1];
rz(-0.46316163) q[3];
sx q[3];
rz(-2.4235631) q[3];
sx q[3];
rz(-2.5754186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51155382) q[2];
sx q[2];
rz(-1.0147107) q[2];
sx q[2];
rz(-2.3021728) q[2];
rz(1.9073585) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(-3.1101826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3473174) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(-0.41912249) q[0];
rz(1.4979111) q[1];
sx q[1];
rz(-1.3587911) q[1];
sx q[1];
rz(0.43325123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.142665) q[0];
sx q[0];
rz(-1.4523399) q[0];
sx q[0];
rz(-1.1112369) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5574607) q[2];
sx q[2];
rz(-0.91382256) q[2];
sx q[2];
rz(1.8767534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2280088) q[1];
sx q[1];
rz(-0.59803793) q[1];
sx q[1];
rz(-1.0526471) q[1];
rz(-pi) q[2];
rz(2.8178704) q[3];
sx q[3];
rz(-0.60327134) q[3];
sx q[3];
rz(-2.3761185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7822632) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(2.9456054) q[2];
rz(2.0187142) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.6897197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48542431) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(-1.0957023) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(-2.9313472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4097262) q[0];
sx q[0];
rz(-1.7667701) q[0];
sx q[0];
rz(2.1541697) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91598454) q[2];
sx q[2];
rz(-1.4794) q[2];
sx q[2];
rz(-0.86751988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2481459) q[1];
sx q[1];
rz(-2.4706744) q[1];
sx q[1];
rz(2.4144961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59981646) q[3];
sx q[3];
rz(-1.9164403) q[3];
sx q[3];
rz(1.4100072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11253396) q[2];
sx q[2];
rz(-1.5052648) q[2];
sx q[2];
rz(-2.0564334) q[2];
rz(-2.8339913) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(2.4881081) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6899684) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(0.65144173) q[1];
sx q[1];
rz(-2.1919498) q[1];
sx q[1];
rz(0.776074) q[1];
rz(-1.5522926) q[2];
sx q[2];
rz(-0.69102364) q[2];
sx q[2];
rz(-0.24204248) q[2];
rz(-0.93919803) q[3];
sx q[3];
rz(-1.152335) q[3];
sx q[3];
rz(-1.4622968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
