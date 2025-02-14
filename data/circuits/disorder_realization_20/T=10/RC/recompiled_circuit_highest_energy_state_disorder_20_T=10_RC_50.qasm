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
rz(-2.1751997) q[0];
sx q[0];
rz(-1.3858495) q[0];
sx q[0];
rz(0.59544271) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(-0.59102494) q[1];
sx q[1];
rz(0.12330595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7656674) q[0];
sx q[0];
rz(-1.0320837) q[0];
sx q[0];
rz(-0.49077351) q[0];
rz(0.78211981) q[2];
sx q[2];
rz(-0.95172666) q[2];
sx q[2];
rz(-1.9591046) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1733381) q[1];
sx q[1];
rz(-2.9072793) q[1];
sx q[1];
rz(-1.1543177) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38148613) q[3];
sx q[3];
rz(-2.7499928) q[3];
sx q[3];
rz(-0.42575437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82723242) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(-0.11191351) q[2];
rz(1.7498451) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258485) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(0.14990212) q[0];
rz(0.28597486) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(1.1289977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47685806) q[0];
sx q[0];
rz(-1.5405415) q[0];
sx q[0];
rz(-1.2824351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.01458077) q[2];
sx q[2];
rz(-1.6304071) q[2];
sx q[2];
rz(-2.4042442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4399229) q[1];
sx q[1];
rz(-1.7365125) q[1];
sx q[1];
rz(-0.36335899) q[1];
x q[2];
rz(-2.1993447) q[3];
sx q[3];
rz(-1.0744922) q[3];
sx q[3];
rz(1.7402349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4631606) q[2];
sx q[2];
rz(-0.9223991) q[2];
sx q[2];
rz(1.9909667) q[2];
rz(-1.1848508) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(-3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48035574) q[0];
sx q[0];
rz(-2.895597) q[0];
sx q[0];
rz(-2.7599957) q[0];
rz(-1.8474139) q[1];
sx q[1];
rz(-2.0353863) q[1];
sx q[1];
rz(2.9433184) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7679032) q[0];
sx q[0];
rz(-1.3438017) q[0];
sx q[0];
rz(-0.17902184) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90132331) q[2];
sx q[2];
rz(-1.633989) q[2];
sx q[2];
rz(2.2128999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4130062) q[1];
sx q[1];
rz(-2.8236514) q[1];
sx q[1];
rz(-2.4012662) q[1];
rz(1.4573077) q[3];
sx q[3];
rz(-0.17593613) q[3];
sx q[3];
rz(2.4324377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3759489) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(-0.51885968) q[2];
rz(1.2505924) q[3];
sx q[3];
rz(-2.0542681) q[3];
sx q[3];
rz(-2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4984703) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(-0.44922391) q[0];
rz(-1.4462224) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(-0.62072388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6889544) q[0];
sx q[0];
rz(-1.7151378) q[0];
sx q[0];
rz(-2.4792433) q[0];
rz(-1.6942107) q[2];
sx q[2];
rz(-2.2006249) q[2];
sx q[2];
rz(0.067867756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6606969) q[1];
sx q[1];
rz(-0.45240739) q[1];
sx q[1];
rz(2.6651732) q[1];
rz(-pi) q[2];
rz(-1.4820547) q[3];
sx q[3];
rz(-0.54026287) q[3];
sx q[3];
rz(1.0080573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(2.980496) q[2];
rz(-2.7437239) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(-2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905269) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.2535198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.391025) q[0];
sx q[0];
rz(-1.5526958) q[0];
sx q[0];
rz(-1.5347634) q[0];
rz(-pi) q[1];
rz(0.80998151) q[2];
sx q[2];
rz(-1.4000386) q[2];
sx q[2];
rz(-2.8851938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67318007) q[1];
sx q[1];
rz(-0.77408964) q[1];
sx q[1];
rz(1.5062529) q[1];
x q[2];
rz(0.0048794172) q[3];
sx q[3];
rz(-2.7944428) q[3];
sx q[3];
rz(2.7018869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72838655) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(1.8750635) q[2];
rz(2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(-0.5411886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7414339) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(3.1165282) q[0];
rz(1.2114245) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-0.2549583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1392649) q[0];
sx q[0];
rz(-3.0920017) q[0];
sx q[0];
rz(0.10164405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0588875) q[2];
sx q[2];
rz(-2.1016663) q[2];
sx q[2];
rz(1.4029897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3820262) q[1];
sx q[1];
rz(-2.3236378) q[1];
sx q[1];
rz(2.3809529) q[1];
rz(1.8942157) q[3];
sx q[3];
rz(-2.5445523) q[3];
sx q[3];
rz(-0.14123329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(-2.6825405) q[2];
rz(-1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727305) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(3.0555308) q[0];
rz(-2.22279) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.930621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5945307) q[0];
sx q[0];
rz(-1.8939202) q[0];
sx q[0];
rz(-1.7842152) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44081273) q[2];
sx q[2];
rz(-2.0939504) q[2];
sx q[2];
rz(0.8280821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43566464) q[1];
sx q[1];
rz(-1.8907232) q[1];
sx q[1];
rz(-0.21448696) q[1];
x q[2];
rz(1.2065229) q[3];
sx q[3];
rz(-1.9134054) q[3];
sx q[3];
rz(-2.2628257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3226037) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(-2.4244579) q[2];
rz(1.0273733) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-2.7023442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913919) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(-0.014658654) q[0];
rz(-0.75153366) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(-1.6709447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5605222) q[0];
sx q[0];
rz(-1.7092686) q[0];
sx q[0];
rz(-0.19000368) q[0];
x q[1];
rz(-2.8445817) q[2];
sx q[2];
rz(-0.26784617) q[2];
sx q[2];
rz(1.0849407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9110221) q[1];
sx q[1];
rz(-1.9729905) q[1];
sx q[1];
rz(-3.0190617) q[1];
rz(-pi) q[2];
rz(-1.5812381) q[3];
sx q[3];
rz(-0.34509788) q[3];
sx q[3];
rz(1.0853678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26059255) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(1.986844) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(0.65034136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(-1.0754841) q[0];
rz(0.12818809) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(0.34559616) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3227279) q[0];
sx q[0];
rz(-0.68536192) q[0];
sx q[0];
rz(-1.8249874) q[0];
rz(-pi) q[1];
rz(1.0900201) q[2];
sx q[2];
rz(-0.62453237) q[2];
sx q[2];
rz(2.9337954) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84050485) q[1];
sx q[1];
rz(-1.6337758) q[1];
sx q[1];
rz(0.54160629) q[1];
rz(-1.4050499) q[3];
sx q[3];
rz(-2.1624544) q[3];
sx q[3];
rz(-1.6935284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(-1.6205988) q[2];
rz(1.3711551) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(0.41485205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82776752) q[0];
sx q[0];
rz(-2.1791552) q[0];
sx q[0];
rz(-2.9058822) q[0];
rz(0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(1.3689573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9873857) q[0];
sx q[0];
rz(-1.1416178) q[0];
sx q[0];
rz(-2.3230419) q[0];
rz(2.1718352) q[2];
sx q[2];
rz(-2.1548993) q[2];
sx q[2];
rz(-0.36512185) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5342218) q[1];
sx q[1];
rz(-1.0080308) q[1];
sx q[1];
rz(0.20197473) q[1];
rz(-2.7810516) q[3];
sx q[3];
rz(-1.4829794) q[3];
sx q[3];
rz(-0.58780625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77091757) q[2];
sx q[2];
rz(-1.8149899) q[2];
sx q[2];
rz(-3.0885922) q[2];
rz(0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5175405) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(-1.3507631) q[1];
sx q[1];
rz(-2.7443934) q[1];
sx q[1];
rz(2.9846356) q[1];
rz(-2.4443632) q[2];
sx q[2];
rz(-1.6710499) q[2];
sx q[2];
rz(1.2895126) q[2];
rz(-2.5637473) q[3];
sx q[3];
rz(-2.2457849) q[3];
sx q[3];
rz(-0.59413322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
