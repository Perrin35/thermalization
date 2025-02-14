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
rz(0.97839657) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7534417) q[0];
sx q[0];
rz(-1.7072258) q[0];
sx q[0];
rz(0.0078972422) q[0];
x q[1];
rz(-0.70887237) q[2];
sx q[2];
rz(-1.1507431) q[2];
sx q[2];
rz(1.5839029) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.96127993) q[1];
sx q[1];
rz(-0.35113564) q[1];
sx q[1];
rz(1.1457972) q[1];
rz(-pi) q[2];
rz(3.0618598) q[3];
sx q[3];
rz(-1.0736205) q[3];
sx q[3];
rz(-0.91547188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3413099) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(2.4260803) q[2];
rz(1.7320775) q[3];
sx q[3];
rz(-2.2655497) q[3];
sx q[3];
rz(-1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187123) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(3.1067644) q[0];
rz(-0.71827978) q[1];
sx q[1];
rz(-1.7363345) q[1];
sx q[1];
rz(-1.9504331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1615526) q[0];
sx q[0];
rz(-1.411952) q[0];
sx q[0];
rz(0.94670403) q[0];
rz(-pi) q[1];
rz(0.18152512) q[2];
sx q[2];
rz(-2.4511538) q[2];
sx q[2];
rz(1.7275563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9893966) q[1];
sx q[1];
rz(-0.65534822) q[1];
sx q[1];
rz(-0.46391497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9510849) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(-1.0781168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8027773) q[2];
sx q[2];
rz(-1.2884527) q[2];
sx q[2];
rz(-3.0954933) q[2];
rz(0.80373803) q[3];
sx q[3];
rz(-1.9019144) q[3];
sx q[3];
rz(2.5900335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(1.4285202) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(-2.5285316) q[0];
rz(2.8763981) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(-0.048096098) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10079028) q[0];
sx q[0];
rz(-1.5682966) q[0];
sx q[0];
rz(1.0889755) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71718054) q[2];
sx q[2];
rz(-2.4240026) q[2];
sx q[2];
rz(0.24914514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7386507) q[1];
sx q[1];
rz(-1.623282) q[1];
sx q[1];
rz(1.5110635) q[1];
rz(-pi) q[2];
rz(-0.35086029) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(-1.3223437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4464104) q[2];
sx q[2];
rz(-1.9789275) q[2];
sx q[2];
rz(-0.41333684) q[2];
rz(2.4731596) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(-1.7641164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9628825) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(2.9982153) q[0];
rz(0.42477056) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(2.7075148) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2370983) q[0];
sx q[0];
rz(-1.6088271) q[0];
sx q[0];
rz(-1.4774986) q[0];
x q[1];
rz(-2.9888568) q[2];
sx q[2];
rz(-2.483791) q[2];
sx q[2];
rz(1.7694103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.078316704) q[1];
sx q[1];
rz(-2.6022291) q[1];
sx q[1];
rz(-1.8453979) q[1];
x q[2];
rz(-0.20695417) q[3];
sx q[3];
rz(-0.82217875) q[3];
sx q[3];
rz(-1.5311629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8726337) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(0.41716519) q[2];
rz(-2.4560691) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(-3.0912257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79614574) q[0];
sx q[0];
rz(-2.8247927) q[0];
sx q[0];
rz(1.6337122) q[0];
rz(2.4729074) q[1];
sx q[1];
rz(-1.1983127) q[1];
sx q[1];
rz(2.0100458) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9848081) q[0];
sx q[0];
rz(-0.97301345) q[0];
sx q[0];
rz(-3.1020107) q[0];
rz(0.3338694) q[2];
sx q[2];
rz(-1.0354932) q[2];
sx q[2];
rz(-1.4620505) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2315438) q[1];
sx q[1];
rz(-1.1041512) q[1];
sx q[1];
rz(-1.779516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82238166) q[3];
sx q[3];
rz(-2.9677792) q[3];
sx q[3];
rz(1.2863152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23182997) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(-0.46169272) q[2];
rz(-0.26990226) q[3];
sx q[3];
rz(-1.88068) q[3];
sx q[3];
rz(-0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5676024) q[0];
sx q[0];
rz(-2.7157768) q[0];
sx q[0];
rz(-1.061576) q[0];
rz(-0.65879446) q[1];
sx q[1];
rz(-1.4219204) q[1];
sx q[1];
rz(-0.40491358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39570582) q[0];
sx q[0];
rz(-0.86942023) q[0];
sx q[0];
rz(1.4331791) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5091841) q[2];
sx q[2];
rz(-0.64130613) q[2];
sx q[2];
rz(-2.8709047) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15952739) q[1];
sx q[1];
rz(-1.7100466) q[1];
sx q[1];
rz(1.0604481) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2085034) q[3];
sx q[3];
rz(-1.8130867) q[3];
sx q[3];
rz(2.1325496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0460661) q[2];
sx q[2];
rz(-1.0045071) q[2];
sx q[2];
rz(-1.7556165) q[2];
rz(-2.4447377) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(0.81010404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5974469) q[0];
sx q[0];
rz(-2.5283041) q[0];
sx q[0];
rz(0.59462732) q[0];
rz(1.1320629) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(-2.6304257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88379556) q[0];
sx q[0];
rz(-0.57475677) q[0];
sx q[0];
rz(3.0416328) q[0];
x q[1];
rz(-2.2726502) q[2];
sx q[2];
rz(-0.41837439) q[2];
sx q[2];
rz(1.6623896) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6936924) q[1];
sx q[1];
rz(-1.3096454) q[1];
sx q[1];
rz(3.0413701) q[1];
rz(-2.4092402) q[3];
sx q[3];
rz(-1.9480289) q[3];
sx q[3];
rz(-0.66483745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1749997) q[2];
sx q[2];
rz(-2.9006557) q[2];
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
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6519046) q[0];
sx q[0];
rz(-2.7629485) q[0];
sx q[0];
rz(-0.65377671) q[0];
rz(0.50103029) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(-2.3562866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0609496) q[0];
sx q[0];
rz(-1.4333581) q[0];
sx q[0];
rz(-2.5713443) q[0];
x q[1];
rz(-2.0205604) q[2];
sx q[2];
rz(-1.0503979) q[2];
sx q[2];
rz(-2.2614417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.874215) q[1];
sx q[1];
rz(-1.7804044) q[1];
sx q[1];
rz(-2.4506635) q[1];
x q[2];
rz(-1.4887962) q[3];
sx q[3];
rz(-1.7991711) q[3];
sx q[3];
rz(0.60199753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.952897) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(2.7719899) q[2];
rz(-0.26436198) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(-3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13067506) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(-1.6226907) q[0];
rz(-1.9021289) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(0.94737238) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6560871) q[0];
sx q[0];
rz(-1.6316905) q[0];
sx q[0];
rz(1.1709309) q[0];
x q[1];
rz(0.16504176) q[2];
sx q[2];
rz(-0.36569917) q[2];
sx q[2];
rz(0.62830454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2202962) q[1];
sx q[1];
rz(-1.4054757) q[1];
sx q[1];
rz(0.16311793) q[1];
rz(-pi) q[2];
rz(-1.9188493) q[3];
sx q[3];
rz(-0.17379856) q[3];
sx q[3];
rz(-0.17979515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47306791) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(-1.1119941) q[2];
rz(0.2019349) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(-0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875459) q[0];
sx q[0];
rz(-2.9534979) q[0];
sx q[0];
rz(1.4659708) q[0];
rz(1.6000043) q[1];
sx q[1];
rz(-1.3009289) q[1];
sx q[1];
rz(-2.8864554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2422243) q[0];
sx q[0];
rz(-1.5535019) q[0];
sx q[0];
rz(-1.7045185) q[0];
x q[1];
rz(2.3681247) q[2];
sx q[2];
rz(-1.8102526) q[2];
sx q[2];
rz(-2.4635893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49225351) q[1];
sx q[1];
rz(-1.6494177) q[1];
sx q[1];
rz(-1.2921974) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19740634) q[3];
sx q[3];
rz(-2.3623086) q[3];
sx q[3];
rz(-3.0912661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5280891) q[2];
sx q[2];
rz(-1.9518423) q[2];
sx q[2];
rz(-0.36821723) q[2];
rz(-0.28299371) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1267283) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(0.49900838) q[1];
sx q[1];
rz(-1.2393163) q[1];
sx q[1];
rz(-2.7543482) q[1];
rz(0.3327315) q[2];
sx q[2];
rz(-1.4495779) q[2];
sx q[2];
rz(-2.4548893) q[2];
rz(-1.3697237) q[3];
sx q[3];
rz(-1.5101505) q[3];
sx q[3];
rz(2.8240614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
