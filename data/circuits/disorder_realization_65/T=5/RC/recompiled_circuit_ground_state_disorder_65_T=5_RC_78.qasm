OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.68706694) q[0];
sx q[0];
rz(-1.878976) q[0];
sx q[0];
rz(2.4071121) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(-1.2999111) q[1];
sx q[1];
rz(-2.1631961) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1837195) q[0];
sx q[0];
rz(-1.5629725) q[0];
sx q[0];
rz(1.70723) q[0];
rz(-pi) q[1];
rz(-0.70887237) q[2];
sx q[2];
rz(-1.1507431) q[2];
sx q[2];
rz(-1.5576897) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96127993) q[1];
sx q[1];
rz(-2.790457) q[1];
sx q[1];
rz(1.1457972) q[1];
x q[2];
rz(3.0618598) q[3];
sx q[3];
rz(-2.0679722) q[3];
sx q[3];
rz(-2.2261208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80028278) q[2];
sx q[2];
rz(-1.6494696) q[2];
sx q[2];
rz(-2.4260803) q[2];
rz(1.4095151) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(-1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8228804) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(-0.034828287) q[0];
rz(0.71827978) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(-1.9504331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6189335) q[0];
sx q[0];
rz(-0.95574035) q[0];
sx q[0];
rz(-2.9466924) q[0];
rz(-1.4227563) q[2];
sx q[2];
rz(-2.2477373) q[2];
sx q[2];
rz(1.1803152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4264448) q[1];
sx q[1];
rz(-2.1472125) q[1];
sx q[1];
rz(1.9020686) q[1];
rz(-pi) q[2];
rz(1.9610522) q[3];
sx q[3];
rz(-2.6724216) q[3];
sx q[3];
rz(0.64599761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8027773) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(-3.0954933) q[2];
rz(2.3378546) q[3];
sx q[3];
rz(-1.2396783) q[3];
sx q[3];
rz(2.5900335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.7130724) q[0];
sx q[0];
rz(-1.5794733) q[0];
sx q[0];
rz(0.61306104) q[0];
rz(-2.8763981) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(-3.0934966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6702794) q[0];
sx q[0];
rz(-1.0889772) q[0];
sx q[0];
rz(-0.0028208931) q[0];
rz(1.0499642) q[2];
sx q[2];
rz(-1.0522868) q[2];
sx q[2];
rz(-1.1074378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16471699) q[1];
sx q[1];
rz(-1.6304468) q[1];
sx q[1];
rz(3.0890134) q[1];
rz(-0.78002413) q[3];
sx q[3];
rz(-0.47545537) q[3];
sx q[3];
rz(-0.96987039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6951822) q[2];
sx q[2];
rz(-1.9789275) q[2];
sx q[2];
rz(2.7282558) q[2];
rz(0.6684331) q[3];
sx q[3];
rz(-1.8139402) q[3];
sx q[3];
rz(-1.7641164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17871019) q[0];
sx q[0];
rz(-0.2186192) q[0];
sx q[0];
rz(-2.9982153) q[0];
rz(2.7168221) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(-2.7075148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804337) q[0];
sx q[0];
rz(-1.4775663) q[0];
sx q[0];
rz(-3.1033959) q[0];
rz(-pi) q[1];
rz(-1.4537895) q[2];
sx q[2];
rz(-2.2196349) q[2];
sx q[2];
rz(-1.9615356) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39555711) q[1];
sx q[1];
rz(-2.0878804) q[1];
sx q[1];
rz(2.9806869) q[1];
x q[2];
rz(-2.9346385) q[3];
sx q[3];
rz(-2.3194139) q[3];
sx q[3];
rz(1.6104298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.268959) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(-0.41716519) q[2];
rz(-2.4560691) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(0.050366966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3454469) q[0];
sx q[0];
rz(-0.31679994) q[0];
sx q[0];
rz(1.6337122) q[0];
rz(-2.4729074) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(-1.1315469) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9145582) q[0];
sx q[0];
rz(-0.59893227) q[0];
sx q[0];
rz(-1.6288487) q[0];
x q[1];
rz(1.0102398) q[2];
sx q[2];
rz(-1.2850637) q[2];
sx q[2];
rz(3.07522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3858429) q[1];
sx q[1];
rz(-1.3846701) q[1];
sx q[1];
rz(2.666074) q[1];
rz(-pi) q[2];
rz(-1.442836) q[3];
sx q[3];
rz(-1.4528414) q[3];
sx q[3];
rz(-1.02533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9097627) q[2];
sx q[2];
rz(-0.54002395) q[2];
sx q[2];
rz(-0.46169272) q[2];
rz(2.8716904) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(-2.4902978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5676024) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(2.0800166) q[0];
rz(-2.4827982) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-0.40491358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.957307) q[0];
sx q[0];
rz(-0.71249639) q[0];
sx q[0];
rz(-2.980583) q[0];
x q[1];
rz(0.63240856) q[2];
sx q[2];
rz(-0.64130613) q[2];
sx q[2];
rz(2.8709047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3337219) q[1];
sx q[1];
rz(-2.0757339) q[1];
sx q[1];
rz(2.9823279) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8431776) q[3];
sx q[3];
rz(-0.95456353) q[3];
sx q[3];
rz(2.7557834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0460661) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(1.7556165) q[2];
rz(-0.69685495) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(2.3314886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5974469) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(-0.59462732) q[0];
rz(-2.0095297) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(-2.6304257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2577971) q[0];
sx q[0];
rz(-0.57475677) q[0];
sx q[0];
rz(-3.0416328) q[0];
rz(-pi) q[1];
rz(-2.862045) q[2];
sx q[2];
rz(-1.8862533) q[2];
sx q[2];
rz(0.73260288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0446544) q[1];
sx q[1];
rz(-1.4739828) q[1];
sx q[1];
rz(-1.3083878) q[1];
rz(-2.6066512) q[3];
sx q[3];
rz(-2.3341093) q[3];
sx q[3];
rz(0.51715467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96659294) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(0.28755507) q[2];
rz(2.0368841) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(0.12565676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6519046) q[0];
sx q[0];
rz(-2.7629485) q[0];
sx q[0];
rz(2.4878159) q[0];
rz(-0.50103029) q[1];
sx q[1];
rz(-2.4153695) q[1];
sx q[1];
rz(0.78530606) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637999) q[0];
sx q[0];
rz(-2.1350088) q[0];
sx q[0];
rz(1.7336506) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56674024) q[2];
sx q[2];
rz(-1.9575685) q[2];
sx q[2];
rz(0.45505986) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.056902) q[1];
sx q[1];
rz(-0.71699079) q[1];
sx q[1];
rz(-2.8193974) q[1];
rz(2.9124746) q[3];
sx q[3];
rz(-1.4909298) q[3];
sx q[3];
rz(0.95019482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1886957) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(-2.7719899) q[2];
rz(0.26436198) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(-0.00082357426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0109176) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(-1.5189019) q[0];
rz(1.2394637) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(0.94737238) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0305875) q[0];
sx q[0];
rz(-1.1717142) q[0];
sx q[0];
rz(-3.0754979) q[0];
rz(1.6336254) q[2];
sx q[2];
rz(-1.9312973) q[2];
sx q[2];
rz(-2.6897827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2202962) q[1];
sx q[1];
rz(-1.7361169) q[1];
sx q[1];
rz(0.16311793) q[1];
x q[2];
rz(-1.7343643) q[3];
sx q[3];
rz(-1.6298098) q[3];
sx q[3];
rz(1.407377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47306791) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(2.0295985) q[2];
rz(-0.2019349) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428381) q[0];
sx q[0];
rz(-0.18809479) q[0];
sx q[0];
rz(1.6756219) q[0];
rz(1.5415883) q[1];
sx q[1];
rz(-1.3009289) q[1];
sx q[1];
rz(-0.25513729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8993683) q[0];
sx q[0];
rz(-1.5535019) q[0];
sx q[0];
rz(1.4370741) q[0];
x q[1];
rz(2.8053984) q[2];
sx q[2];
rz(-2.3393235) q[2];
sx q[2];
rz(-0.65450571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49225351) q[1];
sx q[1];
rz(-1.6494177) q[1];
sx q[1];
rz(1.2921974) q[1];
rz(0.76948036) q[3];
sx q[3];
rz(-1.7090685) q[3];
sx q[3];
rz(-1.4797803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61350358) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(-2.7733754) q[2];
rz(0.28299371) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(-0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.014864347) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(0.49900838) q[1];
sx q[1];
rz(-1.2393163) q[1];
sx q[1];
rz(-2.7543482) q[1];
rz(-1.4426184) q[2];
sx q[2];
rz(-1.9009931) q[2];
sx q[2];
rz(-0.92585678) q[2];
rz(-1.7718689) q[3];
sx q[3];
rz(-1.6314421) q[3];
sx q[3];
rz(-0.31753123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
