OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(0.84375381) q[0];
rz(1.002797) q[1];
sx q[1];
rz(-2.306814) q[1];
sx q[1];
rz(0.53258449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7850936) q[0];
sx q[0];
rz(-0.66996941) q[0];
sx q[0];
rz(1.8125303) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5101846) q[2];
sx q[2];
rz(-1.4483588) q[2];
sx q[2];
rz(0.40529006) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4733359) q[1];
sx q[1];
rz(-1.7219667) q[1];
sx q[1];
rz(-2.2676047) q[1];
rz(-pi) q[2];
rz(0.37181969) q[3];
sx q[3];
rz(-1.9247247) q[3];
sx q[3];
rz(1.8077426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3089932) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(0.12283202) q[2];
rz(3.099856) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(-2.6532555) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9915344) q[0];
sx q[0];
rz(-0.3287065) q[0];
sx q[0];
rz(1.1154037) q[0];
rz(0.54488048) q[1];
sx q[1];
rz(-1.7439525) q[1];
sx q[1];
rz(-2.5741408) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7780964) q[0];
sx q[0];
rz(-1.6045185) q[0];
sx q[0];
rz(3.0959227) q[0];
x q[1];
rz(1.602515) q[2];
sx q[2];
rz(-0.79481541) q[2];
sx q[2];
rz(1.0358126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48914805) q[1];
sx q[1];
rz(-0.57935435) q[1];
sx q[1];
rz(2.9577319) q[1];
rz(-0.48226815) q[3];
sx q[3];
rz(-1.8390788) q[3];
sx q[3];
rz(-2.0577459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51022092) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(0.1046293) q[2];
rz(-2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-0.73205718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0013244) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.6562756) q[0];
rz(-2.5775919) q[1];
sx q[1];
rz(-2.0878849) q[1];
sx q[1];
rz(-2.3174813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38496537) q[0];
sx q[0];
rz(-0.55776309) q[0];
sx q[0];
rz(-1.2560448) q[0];
x q[1];
rz(2.3468909) q[2];
sx q[2];
rz(-1.1247171) q[2];
sx q[2];
rz(0.54903713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4802758) q[1];
sx q[1];
rz(-1.3696028) q[1];
sx q[1];
rz(2.7072705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7830267) q[3];
sx q[3];
rz(-1.5023059) q[3];
sx q[3];
rz(-2.2048782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.04909) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.606344) q[2];
rz(-2.2384079) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(-2.7170392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20151888) q[0];
sx q[0];
rz(-0.98648447) q[0];
sx q[0];
rz(2.7281813) q[0];
rz(-2.9460733) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(1.6874541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.337404) q[0];
sx q[0];
rz(-1.4436663) q[0];
sx q[0];
rz(-2.26343) q[0];
rz(-pi) q[1];
rz(2.7858673) q[2];
sx q[2];
rz(-1.3918096) q[2];
sx q[2];
rz(2.9236874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6100055) q[1];
sx q[1];
rz(-2.0466699) q[1];
sx q[1];
rz(0.31756084) q[1];
rz(-pi) q[2];
rz(0.24202427) q[3];
sx q[3];
rz(-0.70902642) q[3];
sx q[3];
rz(2.5828862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2426408) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(1.2658489) q[2];
rz(2.2856581) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(-1.2148733) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4796546) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(-1.4187752) q[0];
rz(1.7105557) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(0.72351825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531401) q[0];
sx q[0];
rz(-2.5306764) q[0];
sx q[0];
rz(0.31129654) q[0];
x q[1];
rz(-2.3143598) q[2];
sx q[2];
rz(-0.9632573) q[2];
sx q[2];
rz(0.22996685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0192467) q[1];
sx q[1];
rz(-2.2775688) q[1];
sx q[1];
rz(0.0031696894) q[1];
rz(-pi) q[2];
rz(1.5635023) q[3];
sx q[3];
rz(-1.9230584) q[3];
sx q[3];
rz(-1.1895113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7001223) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(2.7777242) q[2];
rz(1.8338592) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(-1.3522805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42809197) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(2.2416903) q[0];
rz(-1.046754) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(0.55211639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22495843) q[0];
sx q[0];
rz(-0.33946589) q[0];
sx q[0];
rz(3.0566264) q[0];
x q[1];
rz(2.5234473) q[2];
sx q[2];
rz(-2.3380438) q[2];
sx q[2];
rz(-1.3230363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99191715) q[1];
sx q[1];
rz(-0.19405716) q[1];
sx q[1];
rz(-2.0582576) q[1];
rz(-0.12142809) q[3];
sx q[3];
rz(-0.59372205) q[3];
sx q[3];
rz(0.93103257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8206574) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(-2.0123539) q[3];
sx q[3];
rz(-1.4720474) q[3];
sx q[3];
rz(-0.78013295) q[3];
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
rz(2.8010913) q[0];
sx q[0];
rz(-1.8393562) q[0];
sx q[0];
rz(0.71402016) q[0];
rz(-1.2227614) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(-0.27539918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3472766) q[0];
sx q[0];
rz(-1.9660608) q[0];
sx q[0];
rz(-2.9121132) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8267639) q[2];
sx q[2];
rz(-0.98199516) q[2];
sx q[2];
rz(2.3662629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95412817) q[1];
sx q[1];
rz(-2.4373332) q[1];
sx q[1];
rz(-2.3440948) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5573165) q[3];
sx q[3];
rz(-2.2960141) q[3];
sx q[3];
rz(2.0303371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068646) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(1.6738221) q[2];
rz(1.5214527) q[3];
sx q[3];
rz(-2.455267) q[3];
sx q[3];
rz(0.93792382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15243212) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(-0.96617019) q[0];
rz(0.78978157) q[1];
sx q[1];
rz(-0.25895324) q[1];
sx q[1];
rz(-2.1678179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6286929) q[0];
sx q[0];
rz(-1.8159465) q[0];
sx q[0];
rz(0.12314491) q[0];
x q[1];
rz(-2.439271) q[2];
sx q[2];
rz(-2.6036545) q[2];
sx q[2];
rz(0.59662102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57339797) q[1];
sx q[1];
rz(-2.1713082) q[1];
sx q[1];
rz(-0.37094122) q[1];
rz(-pi) q[2];
rz(1.4548746) q[3];
sx q[3];
rz(-2.6629857) q[3];
sx q[3];
rz(2.6461864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3070273) q[2];
sx q[2];
rz(-1.7119188) q[2];
sx q[2];
rz(-2.5359421) q[2];
rz(1.0904795) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(-3.0428913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379446) q[0];
sx q[0];
rz(-3.0927959) q[0];
sx q[0];
rz(-0.57156372) q[0];
rz(0.38772186) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-0.8943843) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69686517) q[0];
sx q[0];
rz(-1.0859216) q[0];
sx q[0];
rz(0.98486017) q[0];
x q[1];
rz(-1.3887843) q[2];
sx q[2];
rz(-2.5669206) q[2];
sx q[2];
rz(-2.4383309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91796658) q[1];
sx q[1];
rz(-1.3296179) q[1];
sx q[1];
rz(-1.9832703) q[1];
rz(-0.13988755) q[3];
sx q[3];
rz(-2.6287946) q[3];
sx q[3];
rz(2.1751805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1423219) q[2];
sx q[2];
rz(-0.88185507) q[2];
sx q[2];
rz(0.87232653) q[2];
rz(3.0669323) q[3];
sx q[3];
rz(-1.229769) q[3];
sx q[3];
rz(2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1416624) q[0];
sx q[0];
rz(-2.7003728) q[0];
sx q[0];
rz(2.4040851) q[0];
rz(2.4420338) q[1];
sx q[1];
rz(-1.0195426) q[1];
sx q[1];
rz(-2.2950744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695037) q[0];
sx q[0];
rz(-2.3583625) q[0];
sx q[0];
rz(2.5602719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9969127) q[2];
sx q[2];
rz(-1.5240545) q[2];
sx q[2];
rz(-2.5747204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10900765) q[1];
sx q[1];
rz(-1.4334588) q[1];
sx q[1];
rz(-1.1395543) q[1];
rz(-pi) q[2];
rz(0.40007253) q[3];
sx q[3];
rz(-1.3988675) q[3];
sx q[3];
rz(1.4330869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49125853) q[2];
sx q[2];
rz(-1.9776055) q[2];
sx q[2];
rz(2.7872046) q[2];
rz(-0.77704159) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(-1.0907772) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95018321) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(-1.3632111) q[1];
sx q[1];
rz(-0.53047219) q[1];
sx q[1];
rz(-0.06906876) q[1];
rz(-1.6949881) q[2];
sx q[2];
rz(-0.71659452) q[2];
sx q[2];
rz(2.334395) q[2];
rz(-0.94181413) q[3];
sx q[3];
rz(-0.54472803) q[3];
sx q[3];
rz(1.7036078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
