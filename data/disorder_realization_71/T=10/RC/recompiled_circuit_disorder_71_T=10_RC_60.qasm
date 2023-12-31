OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60183817) q[0];
sx q[0];
rz(-1.9617426) q[0];
sx q[0];
rz(1.3215617) q[0];
x q[1];
rz(-1.285577) q[2];
sx q[2];
rz(-2.5370295) q[2];
sx q[2];
rz(-0.7315469) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0383366) q[1];
sx q[1];
rz(-1.5557319) q[1];
sx q[1];
rz(-1.9008093) q[1];
rz(-pi) q[2];
x q[2];
rz(1.417744) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-2.412964) q[2];
rz(-2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.306863) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(-2.2671525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782198) q[0];
sx q[0];
rz(-1.997943) q[0];
sx q[0];
rz(0.33421974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99526694) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(1.0075943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2494617) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-2.4888121) q[1];
rz(1.7631093) q[3];
sx q[3];
rz(-2.8674539) q[3];
sx q[3];
rz(0.4916693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(1.0748192) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-0.39594617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031035) q[0];
sx q[0];
rz(-1.9651411) q[0];
sx q[0];
rz(2.1179384) q[0];
rz(-pi) q[1];
rz(0.01404889) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(0.77169466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0242651) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(2.4584103) q[1];
rz(-pi) q[2];
rz(2.8860502) q[3];
sx q[3];
rz(-1.5934172) q[3];
sx q[3];
rz(-0.25657755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48196402) q[0];
sx q[0];
rz(-0.11419645) q[0];
sx q[0];
rz(-2.1301079) q[0];
rz(0.31077023) q[2];
sx q[2];
rz(-2.4756873) q[2];
sx q[2];
rz(-3.0096465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.213821) q[1];
sx q[1];
rz(-1.908761) q[1];
sx q[1];
rz(2.4073699) q[1];
x q[2];
rz(2.7531284) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(0.24211611) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2634537) q[0];
sx q[0];
rz(-1.3100776) q[0];
sx q[0];
rz(3.1104452) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3741578) q[2];
sx q[2];
rz(-2.8637297) q[2];
sx q[2];
rz(2.1674736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6619751) q[1];
sx q[1];
rz(-1.2503337) q[1];
sx q[1];
rz(1.5451876) q[1];
rz(-pi) q[2];
rz(1.2916336) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(-2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(-0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(-2.81566) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78544261) q[0];
sx q[0];
rz(-0.66473648) q[0];
sx q[0];
rz(0.63389969) q[0];
x q[1];
rz(-0.58165254) q[2];
sx q[2];
rz(-2.622421) q[2];
sx q[2];
rz(-0.66228629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1833916) q[1];
sx q[1];
rz(-0.32918731) q[1];
sx q[1];
rz(-2.7126461) q[1];
x q[2];
rz(-2.9767354) q[3];
sx q[3];
rz(-3.0391209) q[3];
sx q[3];
rz(-0.41302478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(0.72189271) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(3.022335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0605474) q[0];
sx q[0];
rz(-1.3067055) q[0];
sx q[0];
rz(0.72780769) q[0];
x q[1];
rz(1.4085521) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(2.721399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99265656) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(-2.0723144) q[1];
rz(0.42640949) q[3];
sx q[3];
rz(-2.3873513) q[3];
sx q[3];
rz(2.756556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(1.7822441) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-2.2559821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6016156) q[0];
sx q[0];
rz(-0.5479387) q[0];
sx q[0];
rz(-2.6849296) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1259414) q[2];
sx q[2];
rz(-2.1502697) q[2];
sx q[2];
rz(-1.5660812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13247989) q[1];
sx q[1];
rz(-1.1785893) q[1];
sx q[1];
rz(-2.1177887) q[1];
rz(-pi) q[2];
rz(0.86730154) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-1.0845832) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866859) q[0];
sx q[0];
rz(-2.1011155) q[0];
sx q[0];
rz(-2.9493939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7068107) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(-2.9895363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.218704) q[1];
sx q[1];
rz(-1.522038) q[1];
sx q[1];
rz(-2.7152275) q[1];
x q[2];
rz(2.7157852) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0231126) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4599633) q[0];
sx q[0];
rz(-0.75461331) q[0];
sx q[0];
rz(0.90453903) q[0];
rz(-1.3329266) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(-2.7931917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.69233209) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(-0.87528054) q[1];
x q[2];
rz(2.5246546) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(-2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(-1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(2.1879788) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-0.52195436) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-0.75688731) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
rz(-2.5064777) q[3];
sx q[3];
rz(-2.1182346) q[3];
sx q[3];
rz(2.3253141) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
