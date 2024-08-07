---
layout: page
permalink: /talks/
title: "Talks"
---

<ul class="post-list">
  {%- for post in site.posts -%}
    {%- if post.categories contains "talks" -%}
	<li>
      {%- assign date_format = site.minima.date_format | default: "%b %-d, %Y" -%}
      <span class="post-meta">{{ post.date | date: date_format }}</span>
      <h3>
        <a class="post-link" href="{{ post.url | relative_url }}">
          {{ post.title | escape }}
        </a>
      </h3>
    </li>
    {%- endif -%}
  {%- endfor -%}
</ul>
