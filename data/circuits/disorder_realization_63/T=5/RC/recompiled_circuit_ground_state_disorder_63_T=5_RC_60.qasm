OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(-1.902268) q[0];
sx q[0];
rz(-1.0898606) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(-1.8884594) q[1];
sx q[1];
rz(0.90322948) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.922674) q[0];
sx q[0];
rz(-2.196133) q[0];
sx q[0];
rz(2.3614285) q[0];
rz(-pi) q[1];
rz(-0.014711424) q[2];
sx q[2];
rz(-0.88764578) q[2];
sx q[2];
rz(2.8840182) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.049202327) q[1];
sx q[1];
rz(-0.70194879) q[1];
sx q[1];
rz(2.5939221) q[1];
rz(-pi) q[2];
rz(1.3525659) q[3];
sx q[3];
rz(-2.6317843) q[3];
sx q[3];
rz(2.9657983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(2.7083) q[2];
rz(-0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-2.8210631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1746154) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(-1.2906661) q[0];
rz(0.67944747) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(-0.89250934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3474059) q[0];
sx q[0];
rz(-1.5457602) q[0];
sx q[0];
rz(1.7480609) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4878534) q[2];
sx q[2];
rz(-1.2481261) q[2];
sx q[2];
rz(1.8316837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50607562) q[1];
sx q[1];
rz(-0.92787023) q[1];
sx q[1];
rz(-1.2256114) q[1];
x q[2];
rz(-0.99231798) q[3];
sx q[3];
rz(-2.7035993) q[3];
sx q[3];
rz(2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9925925) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(2.0089669) q[2];
rz(2.4752786) q[3];
sx q[3];
rz(-1.7814813) q[3];
sx q[3];
rz(1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76729524) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(0.98180109) q[0];
rz(0.94217316) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(-2.2312677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8273329) q[0];
sx q[0];
rz(-1.1480486) q[0];
sx q[0];
rz(2.6092922) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40159638) q[2];
sx q[2];
rz(-2.842428) q[2];
sx q[2];
rz(0.49002346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4228607) q[1];
sx q[1];
rz(-1.8721458) q[1];
sx q[1];
rz(1.2217997) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1317838) q[3];
sx q[3];
rz(-1.3512423) q[3];
sx q[3];
rz(-1.9569524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(0.18690898) q[2];
rz(0.16380353) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(-0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(-0.69766587) q[0];
rz(-0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.1950511) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9055525) q[0];
sx q[0];
rz(-2.4396606) q[0];
sx q[0];
rz(-1.0206971) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5279453) q[2];
sx q[2];
rz(-0.23229182) q[2];
sx q[2];
rz(-2.7181546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6704441) q[1];
sx q[1];
rz(-0.19757825) q[1];
sx q[1];
rz(1.0001282) q[1];
rz(-pi) q[2];
rz(3.0680858) q[3];
sx q[3];
rz(-1.2871082) q[3];
sx q[3];
rz(-2.7105041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.4231851) q[2];
sx q[2];
rz(2.9441693) q[2];
rz(-0.10489634) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(-3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5941493) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(1.0092258) q[0];
rz(3.1127473) q[1];
sx q[1];
rz(-1.6639158) q[1];
sx q[1];
rz(0.65863329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26732609) q[0];
sx q[0];
rz(-1.9443041) q[0];
sx q[0];
rz(1.1248571) q[0];
x q[1];
rz(1.9444185) q[2];
sx q[2];
rz(-2.6251617) q[2];
sx q[2];
rz(-0.85631285) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62477905) q[1];
sx q[1];
rz(-0.1923696) q[1];
sx q[1];
rz(2.6386847) q[1];
x q[2];
rz(0.27401383) q[3];
sx q[3];
rz(-2.4276795) q[3];
sx q[3];
rz(1.4307724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7625526) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(0.23400447) q[2];
rz(1.0903357) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91319084) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(-1.2983904) q[0];
rz(0.78650728) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(-1.7826805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3284861) q[0];
sx q[0];
rz(-2.4755602) q[0];
sx q[0];
rz(-0.12909992) q[0];
rz(-pi) q[1];
rz(0.26769079) q[2];
sx q[2];
rz(-1.7867733) q[2];
sx q[2];
rz(2.9962073) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5740093) q[1];
sx q[1];
rz(-1.5720075) q[1];
sx q[1];
rz(3.1399425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2511061) q[3];
sx q[3];
rz(-1.1292919) q[3];
sx q[3];
rz(-0.46605817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(-0.1725014) q[2];
rz(-2.6532459) q[3];
sx q[3];
rz(-1.9988873) q[3];
sx q[3];
rz(-2.02777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(2.8016222) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.2350157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.039073) q[0];
sx q[0];
rz(-2.2082941) q[0];
sx q[0];
rz(-2.0387285) q[0];
rz(2.4900774) q[2];
sx q[2];
rz(-2.5861202) q[2];
sx q[2];
rz(3.0998203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3307747) q[1];
sx q[1];
rz(-2.0502809) q[1];
sx q[1];
rz(-1.564784) q[1];
x q[2];
rz(2.7478479) q[3];
sx q[3];
rz(-1.3662158) q[3];
sx q[3];
rz(2.1285299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(0.55533448) q[2];
rz(-2.9739144) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(-0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.70917201) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(1.1322016) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(2.1999377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96977329) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(-2.9540747) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8563111) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(0.51644737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42313566) q[1];
sx q[1];
rz(-0.65538156) q[1];
sx q[1];
rz(1.7792367) q[1];
rz(-pi) q[2];
rz(-1.7835983) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(1.4546284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7534916) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.5388185) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(-3.0901618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0592773) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(0.96784651) q[0];
rz(-1.2713426) q[1];
sx q[1];
rz(-1.2920734) q[1];
sx q[1];
rz(0.06180067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3508575) q[0];
sx q[0];
rz(-0.684832) q[0];
sx q[0];
rz(2.5032296) q[0];
x q[1];
rz(-1.3702622) q[2];
sx q[2];
rz(-1.8248744) q[2];
sx q[2];
rz(0.24607436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7605412) q[1];
sx q[1];
rz(-0.88812056) q[1];
sx q[1];
rz(-1.6939916) q[1];
rz(-pi) q[2];
rz(-2.8093407) q[3];
sx q[3];
rz(-1.2945172) q[3];
sx q[3];
rz(-2.0656697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(-3.0657213) q[2];
rz(-1.9742981) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9645914) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(0.28468537) q[0];
rz(-0.074706569) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(-1.7237192) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95686326) q[0];
sx q[0];
rz(-1.3780344) q[0];
sx q[0];
rz(-0.75988976) q[0];
rz(-pi) q[1];
rz(-2.9921164) q[2];
sx q[2];
rz(-1.0188531) q[2];
sx q[2];
rz(1.1750482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7602882) q[1];
sx q[1];
rz(-2.3592296) q[1];
sx q[1];
rz(0.60790498) q[1];
rz(-pi) q[2];
x q[2];
rz(3.10793) q[3];
sx q[3];
rz(-0.92516781) q[3];
sx q[3];
rz(2.1098304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(1.2791862) q[2];
rz(2.159436) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(-2.2568683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(-2.1263532) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(-1.954409) q[2];
sx q[2];
rz(-2.8593079) q[2];
sx q[2];
rz(1.5310892) q[2];
rz(0.61839907) q[3];
sx q[3];
rz(-2.3000345) q[3];
sx q[3];
rz(-2.2207501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
