OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(-5.6145515) q[1];
sx q[1];
rz(0.86548391) q[1];
sx q[1];
rz(15.794985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60458175) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(1.0010615) q[0];
rz(-0.76724648) q[2];
sx q[2];
rz(-1.6635206) q[2];
sx q[2];
rz(-2.7887432) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4359856) q[1];
sx q[1];
rz(-1.7090624) q[1];
sx q[1];
rz(-2.0558753) q[1];
x q[2];
rz(0.41214715) q[3];
sx q[3];
rz(-1.5282359) q[3];
sx q[3];
rz(-1.3679078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(-2.7677317) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502991) q[0];
sx q[0];
rz(-2.5851558) q[0];
sx q[0];
rz(-2.2668905) q[0];
x q[1];
rz(1.3090918) q[2];
sx q[2];
rz(-1.7644901) q[2];
sx q[2];
rz(1.9270093) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0034069) q[1];
sx q[1];
rz(-2.1669743) q[1];
sx q[1];
rz(-1.2852933) q[1];
x q[2];
rz(-1.5156636) q[3];
sx q[3];
rz(-1.1109567) q[3];
sx q[3];
rz(0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(-2.9648798) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(0.40476558) q[0];
rz(-1.8602712) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.3084897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256516) q[0];
sx q[0];
rz(-2.0309272) q[0];
sx q[0];
rz(-0.56002754) q[0];
x q[1];
rz(-2.6042125) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(2.1539719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8150427) q[1];
sx q[1];
rz(-0.76154852) q[1];
sx q[1];
rz(1.1990859) q[1];
rz(-0.19034068) q[3];
sx q[3];
rz(-0.44572178) q[3];
sx q[3];
rz(-1.0090855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.219316) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999009) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(-1.0571009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(2.7042424) q[0];
x q[1];
rz(-2.1787203) q[2];
sx q[2];
rz(-1.086237) q[2];
sx q[2];
rz(2.059666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5406487) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(0.01407108) q[1];
x q[2];
rz(2.3841342) q[3];
sx q[3];
rz(-1.7372903) q[3];
sx q[3];
rz(-1.7050626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(-0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7856359) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-2.6729029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5422573) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(-0.50962781) q[0];
x q[1];
rz(-1.7519978) q[2];
sx q[2];
rz(-1.8743519) q[2];
sx q[2];
rz(-2.5163577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3846181) q[1];
sx q[1];
rz(-2.9006835) q[1];
sx q[1];
rz(-0.5730281) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0607871) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(-2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(0.86432499) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7140759) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(2.8053455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7173548) q[0];
sx q[0];
rz(-0.80168085) q[0];
sx q[0];
rz(2.5108811) q[0];
rz(-pi) q[1];
rz(2.3669846) q[2];
sx q[2];
rz(-0.75658549) q[2];
sx q[2];
rz(2.2677383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0834004) q[1];
sx q[1];
rz(-2.5501049) q[1];
sx q[1];
rz(0.085043474) q[1];
x q[2];
rz(1.7429966) q[3];
sx q[3];
rz(-0.88168722) q[3];
sx q[3];
rz(2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(0.0011750778) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0703697) q[0];
sx q[0];
rz(-0.56814146) q[0];
sx q[0];
rz(2.2413261) q[0];
x q[1];
rz(3.1292186) q[2];
sx q[2];
rz(-1.3482598) q[2];
sx q[2];
rz(0.48977938) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.104407) q[1];
sx q[1];
rz(-1.6652602) q[1];
sx q[1];
rz(0.027992804) q[1];
rz(-3.0435211) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(-1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(2.337713) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893387) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(-0.071285204) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55401037) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(0.49288921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56402928) q[1];
sx q[1];
rz(-1.659435) q[1];
sx q[1];
rz(3.1313529) q[1];
rz(3.026268) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(-0.61570864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(2.5381952) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-0.75072748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663561) q[0];
sx q[0];
rz(-2.7285517) q[0];
sx q[0];
rz(-0.77126276) q[0];
rz(-1.6523916) q[2];
sx q[2];
rz(-1.8494693) q[2];
sx q[2];
rz(2.1747053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1474485) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(1.2575498) q[1];
rz(1.5446072) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.3903842) q[0];
rz(2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.4601382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4936192) q[0];
sx q[0];
rz(-0.68035347) q[0];
sx q[0];
rz(-0.86662678) q[0];
rz(-pi) q[1];
rz(-1.0667801) q[2];
sx q[2];
rz(-1.5911615) q[2];
sx q[2];
rz(-2.6838944) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32669386) q[1];
sx q[1];
rz(-1.9683451) q[1];
sx q[1];
rz(1.3709929) q[1];
rz(0.7695997) q[3];
sx q[3];
rz(-1.0255314) q[3];
sx q[3];
rz(1.256497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(-1.4617408) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(2.3537221) q[2];
sx q[2];
rz(-1.315016) q[2];
sx q[2];
rz(-1.6223326) q[2];
rz(2.4685728) q[3];
sx q[3];
rz(-1.2613847) q[3];
sx q[3];
rz(-3.1223084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];