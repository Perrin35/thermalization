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
rz(-1.2004381) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5759597) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(0.54310449) q[0];
rz(-1.0317694) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(0.10345085) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9125036) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(-1.9181262) q[1];
x q[2];
rz(-2.3463763) q[3];
sx q[3];
rz(-1.8024369) q[3];
sx q[3];
rz(-0.29602805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.583064) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(2.6057459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-2.147701) q[0];
sx q[0];
rz(-2.9931086) q[0];
rz(-pi) q[1];
rz(2.4321062) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(-0.62765861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6260813) q[1];
sx q[1];
rz(-1.1647845) q[1];
sx q[1];
rz(1.2938234) q[1];
rz(-pi) q[2];
rz(-1.7137394) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(0.95834857) q[2];
rz(-3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41985837) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(-0.025807468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31608554) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(-0.020629701) q[0];
x q[1];
rz(1.259272) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(2.4740919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7681231) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(-0.5975238) q[1];
rz(-pi) q[2];
rz(1.1879376) q[3];
sx q[3];
rz(-1.6930876) q[3];
sx q[3];
rz(2.0910636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(1.5766778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(0.53612971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6484084) q[0];
sx q[0];
rz(-2.3097561) q[0];
sx q[0];
rz(-2.1302845) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(1.1255217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0587479) q[1];
sx q[1];
rz(-1.1943598) q[1];
sx q[1];
rz(-2.46545) q[1];
rz(1.2767775) q[3];
sx q[3];
rz(-2.2259568) q[3];
sx q[3];
rz(0.40757195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.608425) q[0];
sx q[0];
rz(-2.5072917) q[0];
sx q[0];
rz(-1.2483291) q[0];
rz(-pi) q[1];
rz(2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(-1.2793465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16380285) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(2.8603641) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78955663) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(1.1556926) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(-1.0587143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67260325) q[0];
sx q[0];
rz(-0.4840584) q[0];
sx q[0];
rz(-2.6812535) q[0];
rz(-pi) q[1];
rz(-0.52144737) q[2];
sx q[2];
rz(-1.2505184) q[2];
sx q[2];
rz(2.5495868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9204273) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(1.8799053) q[1];
x q[2];
rz(-0.22535725) q[3];
sx q[3];
rz(-2.4212824) q[3];
sx q[3];
rz(-2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6713312) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-2.2999433) q[1];
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
rz(1.3782578) q[2];
sx q[2];
rz(-2.1888869) q[2];
sx q[2];
rz(0.75380177) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6262498) q[1];
sx q[1];
rz(-1.5222349) q[1];
sx q[1];
rz(-1.6471144) q[1];
rz(0.82960415) q[3];
sx q[3];
rz(-1.4083574) q[3];
sx q[3];
rz(-0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324875) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(1.0743272) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0268961) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(2.7483658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8199181) q[1];
sx q[1];
rz(-2.3707317) q[1];
sx q[1];
rz(2.2714771) q[1];
x q[2];
rz(1.9440218) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(1.9762135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-0.94690698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(-1.253771) q[0];
x q[1];
rz(-3.0445381) q[2];
sx q[2];
rz(-1.6104873) q[2];
sx q[2];
rz(1.5904215) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35518256) q[1];
sx q[1];
rz(-1.9043515) q[1];
sx q[1];
rz(-1.6598318) q[1];
x q[2];
rz(-1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(2.5069619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(-2.4972829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(0.3084348) q[0];
rz(0.23004736) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(2.7237797) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5035142) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(-0.7034941) q[1];
rz(0.82018567) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(-1.5370777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(-0.59990668) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.94639287) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(-1.0944081) q[3];
sx q[3];
rz(-2.4088358) q[3];
sx q[3];
rz(0.81523304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
