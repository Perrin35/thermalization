OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(-2.2154007) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2339708) q[0];
sx q[0];
rz(-1.3682433) q[0];
sx q[0];
rz(0.58637182) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63106545) q[2];
sx q[2];
rz(-1.4853146) q[2];
sx q[2];
rz(0.98110547) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3057184) q[1];
sx q[1];
rz(-0.22138518) q[1];
sx q[1];
rz(2.1165044) q[1];
rz(-pi) q[2];
rz(3.0285809) q[3];
sx q[3];
rz(-1.7859308) q[3];
sx q[3];
rz(2.5290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7131876) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.5361319) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40614906) q[0];
sx q[0];
rz(-0.14649728) q[0];
sx q[0];
rz(-2.2551401) q[0];
rz(-pi) q[1];
rz(2.5218711) q[2];
sx q[2];
rz(-1.1148858) q[2];
sx q[2];
rz(-1.0085307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3000852) q[1];
sx q[1];
rz(-1.2282787) q[1];
sx q[1];
rz(3.1373346) q[1];
rz(-pi) q[2];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.4041956) q[3];
sx q[3];
rz(1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7028971) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(0.3749795) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.1876371) q[1];
sx q[1];
rz(1.8050271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61993116) q[0];
sx q[0];
rz(-1.3494028) q[0];
sx q[0];
rz(-3.0822166) q[0];
rz(1.1281563) q[2];
sx q[2];
rz(-2.1146362) q[2];
sx q[2];
rz(1.3904872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7610443) q[1];
sx q[1];
rz(-0.79954631) q[1];
sx q[1];
rz(1.9995081) q[1];
rz(0.41575723) q[3];
sx q[3];
rz(-2.3299179) q[3];
sx q[3];
rz(-0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0388564) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(2.1170199) q[2];
rz(1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(-1.9112446) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(-1.4039325) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4481903) q[0];
sx q[0];
rz(-1.9172137) q[0];
sx q[0];
rz(-1.3643826) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83547445) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(-0.26963216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2851706) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(1.1986033) q[1];
x q[2];
rz(0.085539354) q[3];
sx q[3];
rz(-2.2604058) q[3];
sx q[3];
rz(2.6181108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4924865) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(-0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3605109) q[0];
sx q[0];
rz(-2.0251944) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(1.5083183) q[1];
sx q[1];
rz(-2.1834178) q[1];
sx q[1];
rz(-0.034084592) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7406697) q[0];
sx q[0];
rz(-1.7310843) q[0];
sx q[0];
rz(0.39808654) q[0];
rz(-0.57406117) q[2];
sx q[2];
rz(-1.7957557) q[2];
sx q[2];
rz(3.0846734) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2195697) q[1];
sx q[1];
rz(-1.8551738) q[1];
sx q[1];
rz(-3.1344096) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94200763) q[3];
sx q[3];
rz(-1.1470456) q[3];
sx q[3];
rz(-2.8429192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(-1.1711228) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532582) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.9357392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57731956) q[0];
sx q[0];
rz(-2.3015129) q[0];
sx q[0];
rz(2.8217836) q[0];
rz(-1.4044936) q[2];
sx q[2];
rz(-0.096509343) q[2];
sx q[2];
rz(-1.6627179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1356493) q[1];
sx q[1];
rz(-1.597153) q[1];
sx q[1];
rz(-0.51075682) q[1];
rz(-pi) q[2];
x q[2];
rz(0.033127012) q[3];
sx q[3];
rz(-1.9702385) q[3];
sx q[3];
rz(1.9314194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(0.051606027) q[2];
rz(-0.40766019) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(2.6962962) q[3];
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
rz(-2.5200941) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(0.53692445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250701) q[0];
sx q[0];
rz(-1.5050097) q[0];
sx q[0];
rz(1.487088) q[0];
rz(-0.14649086) q[2];
sx q[2];
rz(-1.8244201) q[2];
sx q[2];
rz(-0.25073642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.010579022) q[1];
sx q[1];
rz(-2.3691633) q[1];
sx q[1];
rz(1.4355245) q[1];
rz(-1.7241931) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(-0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(-0.34240001) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328338) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(0.38761815) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.7128046) q[1];
sx q[1];
rz(-0.33755916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885376) q[0];
sx q[0];
rz(-2.6153784) q[0];
sx q[0];
rz(0.40286581) q[0];
rz(1.6667716) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(-0.98737398) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6768407) q[1];
sx q[1];
rz(-0.58603474) q[1];
sx q[1];
rz(-1.8094443) q[1];
rz(-pi) q[2];
rz(-2.3861305) q[3];
sx q[3];
rz(-1.6172234) q[3];
sx q[3];
rz(-1.69343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(-1.2672651) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.9545814) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-0.22098456) q[0];
rz(-0.9221319) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(1.0029213) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1598338) q[0];
sx q[0];
rz(-0.28875414) q[0];
sx q[0];
rz(0.50555484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6384764) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(-2.2803277) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79473125) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(-2.5882583) q[1];
rz(1.8995729) q[3];
sx q[3];
rz(-0.71328304) q[3];
sx q[3];
rz(-0.79771358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-0.74206918) q[2];
sx q[2];
rz(-2.9830902) q[2];
rz(0.45378271) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-2.2973072) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-0.89231649) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6979881) q[0];
sx q[0];
rz(-0.83570489) q[0];
sx q[0];
rz(0.23905917) q[0];
rz(-pi) q[1];
rz(-1.6962264) q[2];
sx q[2];
rz(-2.4969366) q[2];
sx q[2];
rz(0.21412011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1103507) q[1];
sx q[1];
rz(-0.24734766) q[1];
sx q[1];
rz(0.27242839) q[1];
rz(0.51384135) q[3];
sx q[3];
rz(-2.6548879) q[3];
sx q[3];
rz(2.0972308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3698547) q[2];
sx q[2];
rz(-0.83738804) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(-1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(0.74444509) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201465) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.66172) q[2];
sx q[2];
rz(-2.8488013) q[2];
sx q[2];
rz(-1.6691096) q[2];
rz(1.5020694) q[3];
sx q[3];
rz(-1.5309661) q[3];
sx q[3];
rz(-1.0504709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
