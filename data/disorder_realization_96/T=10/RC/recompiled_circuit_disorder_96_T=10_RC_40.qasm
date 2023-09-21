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
rz(3.0135305) q[0];
sx q[0];
rz(11.749) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(1.9411545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5759597) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(-2.5984882) q[0];
rz(-pi) q[1];
x q[1];
rz(3.009216) q[2];
sx q[2];
rz(-2.1059603) q[2];
sx q[2];
rz(-1.3995427) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9125036) q[1];
sx q[1];
rz(-0.87516811) q[1];
sx q[1];
rz(1.2234664) q[1];
rz(-pi) q[2];
rz(0.31906268) q[3];
sx q[3];
rz(-0.82108077) q[3];
sx q[3];
rz(2.0882437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(0.53584677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8230096) q[0];
sx q[0];
rz(-0.59360015) q[0];
sx q[0];
rz(-1.3472605) q[0];
x q[1];
rz(-1.8230121) q[2];
sx q[2];
rz(-2.2644342) q[2];
sx q[2];
rz(-2.0351621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9745001) q[1];
sx q[1];
rz(-1.8247316) q[1];
sx q[1];
rz(-0.42029917) q[1];
rz(-pi) q[2];
x q[2];
rz(1.348098) q[3];
sx q[3];
rz(-2.9950812) q[3];
sx q[3];
rz(-0.033586249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-0.95834857) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(0.025807468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2496101) q[0];
sx q[0];
rz(-1.5907856) q[0];
sx q[0];
rz(-1.8206527) q[0];
rz(-pi) q[1];
rz(3.0518968) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(-0.87583625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5379996) q[1];
sx q[1];
rz(-0.94121658) q[1];
sx q[1];
rz(-1.0521207) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13173007) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6648383) q[0];
sx q[0];
rz(-1.1676482) q[0];
sx q[0];
rz(-0.82157764) q[0];
x q[1];
rz(2.9085607) q[2];
sx q[2];
rz(-0.2430025) q[2];
sx q[2];
rz(0.88569966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0828447) q[1];
sx q[1];
rz(-1.1943598) q[1];
sx q[1];
rz(-0.67614268) q[1];
rz(-pi) q[2];
rz(1.8648151) q[3];
sx q[3];
rz(-0.91563581) q[3];
sx q[3];
rz(-2.7340207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-1.0513603) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.608425) q[0];
sx q[0];
rz(-0.63430099) q[0];
sx q[0];
rz(1.8932635) q[0];
x q[1];
rz(2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(1.8622461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4823618) q[1];
sx q[1];
rz(-1.8419957) q[1];
sx q[1];
rz(1.2960474) q[1];
rz(-pi) q[2];
rz(2.99302) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(-0.35219231) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-0.56110704) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832702) q[0];
sx q[0];
rz(-1.1407307) q[0];
sx q[0];
rz(-1.3413315) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5541359) q[2];
sx q[2];
rz(-0.60411462) q[2];
sx q[2];
rz(-0.47746745) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42029542) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(-0.23423127) q[1];
x q[2];
rz(0.22535725) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(0.55707896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.51145) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8182939) q[0];
sx q[0];
rz(-0.57384402) q[0];
sx q[0];
rz(-0.42563514) q[0];
x q[1];
rz(2.5146671) q[2];
sx q[2];
rz(-1.4142087) q[2];
sx q[2];
rz(2.2121034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6517666) q[1];
sx q[1];
rz(-3.0511599) q[1];
sx q[1];
rz(2.138278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3119885) q[3];
sx q[3];
rz(-1.7332352) q[3];
sx q[3];
rz(0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4293489) q[0];
sx q[0];
rz(-2.0040383) q[0];
sx q[0];
rz(-2.594125) q[0];
rz(-pi) q[1];
rz(-1.2049963) q[2];
sx q[2];
rz(-1.6779643) q[2];
sx q[2];
rz(-2.0049713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-1.0091262) q[1];
sx q[1];
rz(2.5820877) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1975708) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(-1.9762135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(0.62966627) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-0.94690698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7787951) q[0];
sx q[0];
rz(-2.6914094) q[0];
sx q[0];
rz(-0.74624004) q[0];
x q[1];
rz(1.6106748) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(0.023488451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8967594) q[1];
sx q[1];
rz(-1.6549126) q[1];
sx q[1];
rz(-2.8068078) q[1];
x q[2];
rz(1.1831207) q[3];
sx q[3];
rz(-0.47260731) q[3];
sx q[3];
rz(1.8567059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.8927195) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-2.4972829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.6258437) q[0];
sx q[0];
rz(-0.17382646) q[0];
x q[1];
rz(-1.8337433) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(0.76542379) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5035142) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(2.4380986) q[1];
rz(-pi) q[2];
rz(2.042291) q[3];
sx q[3];
rz(-1.9715371) q[3];
sx q[3];
rz(-2.5221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(1.7462126) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(1.0944081) q[3];
sx q[3];
rz(-0.73275685) q[3];
sx q[3];
rz(-2.3263596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];